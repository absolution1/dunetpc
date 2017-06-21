#ifndef __MSGLOGGER_H__
#define __MSGLOGGER_H__

#include <ostream>
#include <iostream>

namespace dlardaq
{
  class LogStream
  {
  public:
  LogStream(std::ostream& _out, std::string& _msg) : out(_out), msg(_msg){;}
    template<typename T>
      std::ostream& operator<<(const T& v)
      {out<<msg<<v; return out; }
  protected:
    std::ostream& out;
    std::string   msg;
  };
  
  static std::string _strinfo  = "[ \033[22;32mINFO\033[0m  ] ";
  static std::string _strwarn  = "[ \033[01;35mWARN\033[0m  ] ";
  static std::string _strerror = "[ \033[22;31mERROR\033[0m ] ";
  
  static LogStream msg_info(std::cout, _strinfo);
  static LogStream msg_warn(std::cout, _strwarn);
  static LogStream msg_err(std::cerr, _strerror);
}

#endif
