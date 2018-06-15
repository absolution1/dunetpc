#ifndef ALGPARTS_H
#define ALGPARTS_H

#include <vector>

std::vector<short> frugal_pedestal_sigkill(const std::vector<short>& raw_in,
                                           const int lookahead, const int threshold)
{
    short median=raw_in[0];
    std::vector<short> ped(raw_in.size(), 0);
    // Are we updating the median (because we're not within a hit)?
    bool updating=true;
    for(size_t i=0; i<raw_in.size()-lookahead; ++i){
        short s=raw_in[i]; // The current sample
        short sig_cand=raw_in[i+lookahead]; // The sample a little way ahead
        // Do we later go over the threshold?
        bool cand_above=(sig_cand>median+threshold);
        // Are we currently below the threshold?
        bool current_below=(s<median);
        // Is there an upcoming transition over threshold? If so, freeze the pedestal...
        if(updating && cand_above){
            updating=false;
        }
        // ...until we fall below the pedestal again
        if( (!updating) && (current_below)){
            updating=true;
        }
        // Do the frugal streaming if we're not in a hit
        if(updating){
            if(s>median) ++median;
            if(s<median) --median;
        }
        ped[i]=median;
    }
    // Get the last few samples, which we couldn't do in the main loop
    // because we'd go out-of-bounds
    for(size_t i=raw_in.size()-lookahead; i<raw_in.size(); ++i){
        ped[i]=median;
    }
    return ped;
}

std::vector<short> apply_fir_filter(const std::vector<short>& input,
                                    const size_t ntaps, const short* taps)
{
    std::vector<short> filtered(input.size(), 0);
    for(size_t i=0; i<input.size(); ++i){
        for(size_t j=0; j<ntaps; ++j){
            const size_t index=i>j ? i-j : 0;
            filtered[i]+=input[index]*taps[j];
        }
    }
    return filtered;
}



#endif
