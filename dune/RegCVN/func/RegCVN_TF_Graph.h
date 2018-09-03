////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       RegCVNGraph
//// Authors:     Ilsoo Seong - iseong@uci.edu
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
///               P.Plonski,                      from DUNE, WUT, Sept. 2017
////
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
//// modified from tf_graph.h
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef RegCVNGraph_h
#define RegCVNGraph_h

#include <memory>
#include <vector>
#include <string>

namespace tensorflow
{
    class Session;
    class Tensor;
}

namespace tf
{

class RegCVNGraph
{
public:
    static std::unique_ptr<RegCVNGraph> create(const char* graph_file_name, const unsigned int &ninputs, const std::vector<std::string> & outputs = {})
    {
        bool success;
        std::unique_ptr<RegCVNGraph> ptr(new RegCVNGraph(graph_file_name, ninputs, outputs,  success));
        if (success) { return ptr; }
        else { return nullptr; }
    }

    ~RegCVNGraph();

    std::vector<float> run(const std::vector< std::vector<float> > & x);

    // process vector of 3D inputs, return vector of 1D outputs; use all inputs
    // if samples = -1, or only the specified number of first samples
    std::vector< std::vector<float> > run(
	const std::vector< std::vector< std::vector< std::vector<float> > > > & x,
	long long int samples = -1);
    std::vector< std::vector<float> > run(
	const std::vector< std::vector< std::vector< std::vector<float> > > > & x,
        const unsigned int& ninputs,
	long long int samples = -1);

    std::vector< std::vector< float > > run(const tensorflow::Tensor & x);
    std::vector< std::vector< float > > run(const std::vector< tensorflow::Tensor >& x);

private:
    /// Not-throwing constructor.
    RegCVNGraph(const char* graph_file_name, const unsigned int& ninputs, const std::vector<std::string> & outputs, bool & success);

    tensorflow::Session* fSession;
    std::vector<std::string> fInputNames;
    std::vector< std::string > fOutputNames;
};

} // namespace tf

#endif
