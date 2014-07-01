#include <Bpp/App/ApplicationTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/NexusIoTree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include <gsl/gsl_randist.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

using bpp::Node;
using bpp::TreeTemplate;
using bpp::TreeTemplateTools;

/// From http://stackoverflow.com/questions/11073932/dont-print-trailing-delimiter-stream-iterator-c
template <class C>
auto
print(std::ostream& os, const C& c,
      const std::string& delim = std::string(", "),
      const std::string& open_brace = std::string("{"),
      const std::string& close_brace = std::string("}")
     ) -> decltype(std::begin(c), std::end(c), os)
{
    os << open_brace;
    auto i = std::begin(c);
    auto e = std::end(c);
    if (i != e)
    {
        os << *i;
        for (++i; i != e; ++i)
            os << delim << *i;
    }
    os << close_brace;
    return os;
}


int run_main(int argc, char**argv)
{
    bpp::BppApplication empirical_posterior(argc, argv, "empirical-posterior");
    std::map<string, string> params = empirical_posterior.getParams();

    unique_ptr<bpp::Alphabet> alphabet(bpp::SequenceApplicationTools::getAlphabet(params, "", false));

    // Sites
    unique_ptr<bpp::VectorSiteContainer> sites(bpp::SequenceApplicationTools::getSiteContainer(alphabet.get(), params));
    sites.reset(bpp::SequenceApplicationTools::getSitesToAnalyse(*sites, params, "", true, false));
    bpp::SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // Model, rates
    unique_ptr<bpp::SubstitutionModel> model(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabet.get(), sites.get(), params));
    unique_ptr<bpp::DiscreteDistribution> rate_dist(bpp::PhylogeneticsApplicationTools::getRateDistribution(params));

    // Tree
    unique_ptr<bpp::Tree> tree(bpp::PhylogeneticsApplicationTools::getTree(params));

    // Parameter of exponential prior
    double exp_mean = bpp::ApplicationTools::getDoubleParameter("empirical-posterior.exp_mean", params, 0.1, "", false);

    // Number of points
    size_t steps = static_cast<size_t>(bpp::ApplicationTools::getIntParameter("empirical-posterior.steps", params, 1000, "", false));

    // Output
    std::string output_path = bpp::ApplicationTools::getAFilePath("empirical-posterior.output_path", params, true, false);
    std::string new_taxon = bpp::ApplicationTools::getStringParameter("empirical-posterior.new_taxon", params, "", "", true, false);

    ofstream output_fp(output_path);

    output_fp << "branch_length,prior,likelihood,posterior\n";

    const int node_id = tree->getLeafId(new_taxon);
    double min_bl = 1e-6, max_bl = 1.0;
    double step = max_bl / steps;

    for(size_t i = 0; i < steps; i++) {
        const double branch_length = min_bl + i * step;
        tree->setDistanceToFather(node_id, branch_length);
        bpp::RHomogeneousTreeLikelihood orig_like(*tree, *sites, model.get(), rate_dist.get(), true, false);
        orig_like.initialize();
        const double ll = orig_like.getLogLikelihood();
        const double prior = std::log(gsl_ran_exponential_pdf(branch_length, exp_mean));
        output_fp << branch_length << ',' << prior << ',' << ll << ',' << prior + ll << '\n';
    }

    return 0;
}

int main(int argc, char** argv)
{
    try {
        run_main(argc, argv);
    } catch(bpp::Exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
}
