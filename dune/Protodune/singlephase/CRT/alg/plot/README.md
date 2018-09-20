plot Directory

Contains algorithms for formatting plots in (module, channel) pairs differently.  Some of these 
might be trivial forwarding to TH2D, others might not use TH2D at all.  Every algorithm here derives 
from CRT::ChannelView and exposes its' interface.  Intentionally avoiding class templates so that 
downstream (in a software sense) algorithms can use these implementation classes with my Factory 
pattern a la 3DSTNeutrons.  
