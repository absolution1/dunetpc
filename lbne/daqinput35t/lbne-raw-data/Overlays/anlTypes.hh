#ifndef __ANLTYPES_H__
#define __ANLTYPES_H__

//==============================================================================
// Constants
//==============================================================================

#define MAX_BINS			1024	// Maximum bins in histogram
#define MAX_CHANNELS		12		// Maximum number of channels in digitizer
#define MAX_DISPLAY_DATA	256		// Maximum size of LabWindows table
#define MAX_EVENT_DATA		2046	// Maximum size of waveform in an event
#define MAX_EVENTS			50000	// Maximum events to store from digitizer
#define MAX_PATHNAME_LEN	260		// Maximum length of path including null byte
#define MAX_TIMING_MARKS	32		// Maximum number of timing marks to track per event

#define MAX_CTRL_DATA		256
#define MAX_DEVICES			16
#define MAX_LENGTH_NAME		16
#define MAX_LENGTH_DESC		64

namespace SSPDAQ{

//==============================================================================
// Enumerated Constants
//==============================================================================

enum errorConstants {
	errorNoError			= 0,

	// Device Communication Errors
	errorCommConnect		= 1,
	errorCommDisconnect		= 2,
	errorCommDiscover		= 3,
	errorCommReceive		= 4,
	errorCommSend			= 5,
	errorCommReceiveZero	= 6,
	errorCommReceiveTimeout	= 7,
	errorCommSendZero		= 8,
	errorCommSendTimeout	= 9,
	errorCommType			= 10,
	errorCommPurge			= 11,
	errorCommQueue			= 12,
	
	// Device Data Errors
	errorDataConnect		= 101,
	errorDataLength			= 102,
	errorDataPurge			= 103,
	errorDataQueue			= 104,
	errorDataReceive		= 105,
	errorDataTimeout		= 106,
	
	// LBNE Errors
	errorEventTooLarge		= 201,
	errorEventTooSmall		= 202,
	errorEventTooMany		= 203,
	errorEventHeader		= 204
};

enum retryConstants {
	retryOff	= 0,
	retryOn		= 1
};

enum commConstants {
	commNone	= 0,
	commUDP		= 1,
	commUSB		= 2
};

enum commandConstants {
	cmdNone			= 0,
	// Basic Commands
	cmdRead			= 1,
	cmdReadMask		= 2,
	cmdWrite		= 3,
	cmdWriteMask	= 4,
	// Array Commands
	cmdArrayRead	= 5,
	cmdArrayWrite	= 6,
	// Fifo Commands
	cmdFifoRead		= 7,
	cmdFifoWrite	= 8,
	numCommands
};

enum statusConstants {
	statusNoError		= 0,
	statusSendError		= 1,
	statusReceiveError	= 2,
	statusTimeoutError	= 3,
	statusAddressError	= 4,
	statusAlignError	= 5,
	statusCommandError	= 6,
	statusSizeError		= 7,
	statusWriteError	= 8		// Returned if read-only address is written
};

//==============================================================================
// Types
//==============================================================================

struct EventFlags{
	unsigned char	pileup			: 1;
	unsigned char	polarity		: 1;
	unsigned char	offsetReadout	: 1;
	unsigned char	cfdValid		: 1;
	unsigned char	reserved		: 4;
};

union EventStatus{
	unsigned char		flags;
	EventFlags	flag;
};

struct EventHeader {	// NOTE: Group fields are listed from MSB to LSB
	unsigned int	header;				// 0xAAAAAAAA
	unsigned short	length;				// Packet Length in unsigned ints (including header)
	unsigned short	group1;				// Trigger Type, Status Flags, Header Type
	unsigned short	triggerID;			// Trigger ID
	unsigned short	group2;				// Module ID, Channel ID
	unsigned short	timestamp[4];		// External Timestamp
								// Words 0-1 = Clocks since last sync pulse
								// Words 2-3 = Sync pulse count
	unsigned short	peakSumLow;			// Lower 16 bits of Peak Sum
	unsigned short	group3;				// Offset of Peak, Higher 8 bits of Peak Sum
	unsigned short	preriseLow;			// Lower 16 bits of Prerise
	unsigned short	group4;				// Lower 8 bits of integratedSum, Higher 8 bits of Prerise
	unsigned short	intSumHigh;			// Upper 16 bits of integratedSum
	unsigned short	baseline;			// Baseline
	unsigned short	cfdPoint[4];		// CFD Timestamp Interpolation Points
	unsigned short	intTimestamp[4];	// Internal Timestamp
								// Word 0 = Reserved for interpolation
								// Words 1-2 = 48 bit Timestamp
};

struct Event {
	unsigned int		header;			// 0xAAAAAAAA
	unsigned short		packetLength;
	unsigned short		triggerType;
	EventStatus	status;
	unsigned char		headerType;
	unsigned short		triggerID;
	unsigned short		moduleID;
	unsigned short		channelID;
	unsigned int		syncDelay;
	unsigned int		syncCount;
	int			peakSum;
	char		peakTime;
	unsigned int		prerise;
	unsigned int		integratedSum;
	unsigned short		baseline;
	short		cfdPoint[4];
	unsigned short		intTimestamp[4];
	unsigned short		waveformWords;
	unsigned short		waveform[MAX_EVENT_DATA];
};

struct EventInfo {
	int max;
	int min;
	int timingMark[MAX_TIMING_MARKS];
	int timingMarkColor[MAX_TIMING_MARKS];
	double intEnergy;
	double intEnergyB;
};

struct LBNEVars {
	unsigned int	commType;
	unsigned int	connected;
	unsigned int	useNoComm;
	unsigned int	useUSBComm;
	unsigned int	useUDPComm;
	unsigned int	numDevices;
	int		currDevice;
	int		panelMain;
	int		panelDigitizer;
	int		panelGeneric;
	int		panelTest;
};

struct CtrlHeader {
	unsigned int address;
	unsigned int command;
	unsigned int size;
	unsigned int status;
};

struct CtrlPacket {
	CtrlHeader	header;
	unsigned int		data[MAX_CTRL_DATA];
};

}//namespace SSPDAQ
#endif
