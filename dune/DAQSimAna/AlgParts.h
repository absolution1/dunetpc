#ifndef ALGPARTS_H
#define ALGPARTS_H

#include <vector>

std::vector<short> frugal_pedestal_sigkill(const std::vector<short>& raw_in,
                                           const int lookahead, const int threshold,
                                           const int ncontig)
{
    short median=raw_in[0];
    std::vector<short> ped(raw_in.size(), 0);
    // Are we updating the median (because we're not within a hit)?
    bool updating=true;
    for(size_t i=0; i<raw_in.size()-(lookahead+ncontig); ++i){
        short s=raw_in[i]; // The current sample

        // Look at `ncontig` samples, starting `lookahead` ticks ahead. Are they all above threshold?
        bool all_cands_above=true;
        for(size_t icand=i+lookahead; icand<i+lookahead+ncontig; ++icand){
            if(raw_in[icand]<median+threshold){
                all_cands_above=false;
                break;
            }
        }
        // Are we currently below the threshold?
        bool current_below=(s<median+threshold);
        // Is there an upcoming transition over threshold? If so, freeze the pedestal...
        if(updating && all_cands_above){
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
    for(size_t i=raw_in.size()-(lookahead+ncontig); i<raw_in.size(); ++i){
        ped[i]=median;
    }
    return ped;
}

std::vector<short> frugal_pedestal(const std::vector<short>& raw_in,
                                   const int ncontig)
{
    short median=raw_in[0];
    std::vector<short> ped(raw_in.size(), 0);

    int runningDiff=0;
    for(size_t i=0; i<raw_in.size(); ++i){
        short s=raw_in[i]; // The current sample

        if(s>median) ++runningDiff;
        if(s<median) --runningDiff;

        if(runningDiff > ncontig){
            ++median;
            runningDiff=0;
        }
        if(runningDiff < -1*ncontig){
            --median;
            runningDiff=0;
        }

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
