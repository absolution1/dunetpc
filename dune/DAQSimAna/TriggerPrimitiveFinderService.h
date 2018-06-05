#ifndef TriggerPrimitiveFinderService_h
#define TriggerPrimitiveFinderService_h

#include <vector>
#include <iostream>

class TriggerPrimitiveFinderService {
 
public:
    struct Hit
    {
        Hit(int _channel, int _startTime, int _charge, int _timeOverThreshold)
            : channel(_channel),
              startTime(_startTime),
              charge(_charge),
              timeOverThreshold(_timeOverThreshold)
            {}
        int channel;
        int startTime;
        int charge;
        int timeOverThreshold;
    };

    virtual ~TriggerPrimitiveFinderService() =default;

    virtual std::vector<TriggerPrimitiveFinderService::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples) = 0;
 
};

#ifndef __CLING__
DECLARE_ART_SERVICE_INTERFACE(TriggerPrimitiveFinderService, LEGACY)
#endif

#endif
