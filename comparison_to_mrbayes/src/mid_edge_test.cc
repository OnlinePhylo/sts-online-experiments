#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

//#include <gsl/gsl_randist.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <set>

using namespace std;

using bpp::Node;
using bpp::TreeTemplate;
using bpp::TreeTemplateTools;

vector<unique_ptr<TreeTemplate<Node>>> load_trees(std::map<string, string>& params, string name="")
{
    vector<bpp::Tree*> trees = bpp::PhylogeneticsApplicationTools::getTrees(params, name);
    vector<unique_ptr<TreeTemplate<Node>>> result;
    result.reserve(trees.size());

    for(const bpp::Tree* t : trees) {
        TreeTemplate<Node>* tt = new TreeTemplate<Node>(*t);
        delete t;
        tt->newOutGroup(tt->getLeaves()[0]);
        tt->resetNodesId();
        // Set up root node like we like to
        Node* root = tt->getRootNode();
        const double d = root->getSon(1)->getDistanceToFather();
        root->getSon(0)->setDistanceToFather(root->getSon(0)->getDistanceToFather() + d);
        root->getSon(1)->setDistanceToFather(0.0);
        result.emplace_back(tt);
    }

    return result;
}

vector<Node*> siblings(Node* node) {
    vector<Node*> result;
    if(!node->hasFather())
        return result;

    for(bpp::Node* n : node->getFather()->getSons()) {
        if(n != node)
            result.push_back(n);
    }
    return result;
}

struct Row
{
    size_t tree_id;
    int node_id;
    bool same_as_tree;
    double bl;
    double mid_edge_log_like;
    double mid_edge_log_like_half;
    size_t rank, rank_half;

    string same() const { return same_as_tree ? "yes" : "no"; };
};

int run_main(int argc, char**argv)
{
    bpp::BppApplication compare_topologies(argc, argv, "compare-topologies");
    std::map<string, string> params = compare_topologies.getParams();

    unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));

    // Sites
    unique_ptr<bpp::VectorSiteContainer> sites(bpp::SequenceApplicationTools::getSiteContainer(alphabet.get(), params));
    sites.reset(bpp::SequenceApplicationTools::getSitesToAnalyse(*sites, params, "", true, false));
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Model, rates
    unique_ptr<bpp::SubstitutionModel> model(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), sites.get(), params));
    unique_ptr<bpp::DiscreteDistribution> rate_dist(bpp::PhylogeneticsApplicationTools::getRateDistribution(params));

    // Trees
    vector<unique_ptr<TreeTemplate<Node>>> trees = load_trees(params);

    // Output
    string output_path = bpp::ApplicationTools::getAFilePath("output.file", params, true, false);
    string pruned_taxon = bpp::ApplicationTools::getStringParameter("pruned.taxon", params, "t1");
    ofstream output_file(output_path);
    output_file << "tree,node,same_as_tree,d,mid_edge_log_like,rank,mid_edge_log_like_half,rank_half\n";

    bool quiet = bpp::ApplicationTools::getBooleanParameter("quiet", params, false);
    size_t tree_idx = 0;

    for(unique_ptr<TreeTemplate<Node>>& tree : trees) {
        const int id = tree->getLeafId(pruned_taxon);
        const int father_id = tree->getNode(id)->getFather()->getId();
        const int actual_edge_id = siblings(tree->getNode(id))[0]->getId();
        std::vector<Row> rows;

        // Drop leaf
        TreeTemplateTools::dropLeaf(*tree, pruned_taxon);

        // Add to middle of every edge
        for(const int node_id : tree->getNodesId()) {
            if(tree->getNode(node_id) == tree->getRootNode() || tree->getNode(node_id) == tree->getRootNode()->getSon(1))
                continue;
            // Copy
            TreeTemplate<Node> t(*tree);
            Node *n = t.getNode(node_id);
            Node *father = n->getFather();
            Node *new_node = new Node(father_id);
            Node *new_leaf = new Node(id, pruned_taxon);

            // Update topology
            new_node->addSon(new_leaf);
            father->removeSon(n);
            new_node->addSon(n);
            father->addSon(new_node);

            // Branch lengths
            const double d = n->getDistanceToFather();
            const double ha_d = d / 2;
            new_leaf->setDistanceToFather(0.0);
            new_node->setDistanceToFather(ha_d);
            n->setDistanceToFather(ha_d);

            // Likelihood
            bpp::RHomogeneousTreeLikelihood like(t, *sites, model.get(), rate_dist.get(), true, false);
            like.initialize();
            const double zero_bl_like = like.getLogLikelihood();
            new_leaf->setDistanceToFather(0.5);
            bpp::RHomogeneousTreeLikelihood like2(t, *sites, model.get(), rate_dist.get(), true, false);
            like2.initialize();
            const double half_bl_like = like2.getLogLikelihood();

            rows.push_back(Row{tree_idx, node_id, node_id==actual_edge_id, d, zero_bl_like, half_bl_like, 0, 0});
        }
        std::sort(rows.begin(), rows.end(), [](const Row& r1, const Row& r2) { return r1.mid_edge_log_like > r2.mid_edge_log_like; });
        for(size_t i = 0; i < rows.size(); i++)
            rows[i].rank = i;
        std::sort(rows.begin(), rows.end(), [](const Row& r1, const Row& r2) { return r1.mid_edge_log_like_half > r2.mid_edge_log_like_half; });
        for(size_t i = 0; i < rows.size(); i++)
            rows[i].rank_half = i;
        std::sort(rows.begin(), rows.end(), [](const Row& r1, const Row& r2) { return r1.node_id < r2.node_id ; });
        for(const Row& row : rows) {
            output_file << row.tree_id << ','
                        << row.node_id << ','
                        << row.same() << ','
                        << row.bl << ','
                        << row.mid_edge_log_like << ','
                        << row.rank << ','
                        << row.mid_edge_log_like_half << ','
                        << row.rank_half << '\n';
        }

        tree_idx++;
    }

    return 0;
}

int main(int argc, char** argv)
{
    run_main(argc, argv);
    try {
    } catch(bpp::Exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
}
