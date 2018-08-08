////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       Graph
//// Authors:     Ilsoo Seong - iseong@uci.edu
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
///               P.Plonski,                      from DUNE, WUT, Sept. 2017
////
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
//// modified from tf_graph.h
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Graph_h
#define Graph_h

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

class Graph
{
public:
    static std::unique_ptr<Graph> create(const char* graph_file_name, const unsigned int &ninputs, const std::vector<std::string> & outputs = {})
    {
        bool success;
        std::unique_ptr<Graph> ptr(new Graph(graph_file_name, ninputs, outputs,  success));
        if (success) { return ptr; }
        else { return nullptr; }
    }

    ~Graph();

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
    Graph(const char* graph_file_name, const unsigned int& ninputs, const std::vector<std::string> & outputs, bool & success);

    tensorflow::Session* fSession;
    std::vector<std::string> fInputNames;
    std::vector< std::string > fOutputNames;
};

} // namespace tf

#endif
