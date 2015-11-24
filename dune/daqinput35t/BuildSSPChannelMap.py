#!/usr/bin/python


# Define the (current) hardware mapping
SSP = {}

# Two special cases
SSP[1] = [ (0, 8), (2, 2) ] # CSU + LSU Module
SSP[2] = [ (1, 12) ]        # IU Module

# IU/LBNL modules where SSP=PD Number
for i in [ 3, 5, 7 ]:
    SSP[i] = [ (i, 12) ]

# CSU modules where SSP=PD Number
for i in [ 4, 6 ]:
    SSP[i] = [ (i, 8) ]


# Convert OpDet, Hardware Channel to offline channel
def OfflineChannel(OpDet, HardwareChannel):
    return 12*OpDet + HardwareChannel


# Convert SSP, Hardware Channel to daqheader->group2
def OnlineHeader(ModuleNo, ChannelNo):
    return 16*ModuleNo + ChannelNo

       
for SSPNum in sorted(SSP.keys()):
    SSPChannel = 0
    for OpDet, NChannels in SSP[SSPNum]:
        for PDChannel in range(NChannels):
            offline = OfflineChannel(OpDet, PDChannel)
            online  = OnlineHeader(SSPNum, SSPChannel)
            SSPChannel += 1
            print online, offline

    



