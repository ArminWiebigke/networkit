#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <mutex>

#include "StringBuilder.h"

#ifdef _MSC_VER
	#define NETWORKT_PRETTY_FUNCTION __FUNCSIG__
#else
	#define NETWORKT_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif

#define LOG_LEVEL_FATAL 0
#define LOG_LEVEL_ERROR 1
#define LOG_LEVEL_WARN  2
#define LOG_LEVEL_INFO  3
#define LOG_LEVEL_DEBUG 4
#define LOG_LEVEL_TRACE 5

// For compatibility with existing code, we define the LOG_LEVEL macro.
#if defined(NETWORKIT_RELEASE_LOGGING)
#	define LOG_LEVEL LOG_LEVEL_NONE
#else
#	define LOG_LEVEL LOG_LEVEL_TRACE
#endif

#define FATAL(...) ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::fatal, __VA_ARGS__)
#define FATALF(...) ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::fatal, __VA_ARGS__)

#define ERROR(...) ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::error, __VA_ARGS__)
#define ERRORF(...) ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::error, __VA_ARGS__)

#define WARN(...)  ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::warn,  __VA_ARGS__)
#define WARNF(...)  ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::warn,  __VA_ARGS__)

#define INFO(...)  ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::info,  __VA_ARGS__)
#define INFOF(...)  ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
		::Aux::Log::LogLevel::info,  __VA_ARGS__)

// DEBUG and TRACE are no-ops if NETWORKIT_RELEASE_LOGGING is defined.
#if defined(NETWORKIT_RELEASE_LOGGING)
#	define DEBUG(...) do {} while(false)
#	define DEBUGF(...) do {} while(false)

#	define TRACE(...) do {} while(false)
#	define TRACEF(...) do {} while(false)
#	define TRACEPOINT do {} while(false)
#else
#	define DEBUG(...) ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
			::Aux::Log::LogLevel::debug, __VA_ARGS__)
#	define DEBUGF(...) ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
			::Aux::Log::LogLevel::debug, __VA_ARGS__)

#	define TRACE(...) ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
			::Aux::Log::LogLevel::trace, __VA_ARGS__)
#	define TRACEF(...) ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
			::Aux::Log::LogLevel::trace, __VA_ARGS__)
#	define TRACEPOINT ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__},\
			::Aux::Log::LogLevel::trace, "tracepoint")
#endif // defined(NETWORKIT_RELEASE_LOGGING)

namespace Aux { namespace Log {

struct Location {
	const char* file;
	const char* function;
	const int line;
};

enum class LogLevel {
	trace,
	debug,
	info,
	warn,
	error,
	fatal
};

/**
 * Accept loglevel as string and set.
 * @param logLevel as string
 */
void setLogLevel(std::string logLevel);

/**
 * @return current loglevel as string
 */
std::string getLogLevel();

namespace Settings {

LogLevel getLogLevel();
void setLogLevel(LogLevel p = LogLevel::info);

void setPrintTime(bool b);
bool getPrintTime();

void setPrintLocation(bool b);
bool getPrintLocation();

void setLogfile(const std::string& filename);
}

namespace Impl {
void log(const Location& loc, LogLevel p, const std::string msg);
} //namespace impl

template<typename...T>
void log(const Location& loc, LogLevel p, const T&...args) {
	if(p >= Settings::getLogLevel()) {
		std::stringstream stream;
		printToStream(stream, args...);
		Impl::log(loc, p, stream.str());
	}
}

template<typename...T>
void logF(const Location& loc, LogLevel p, const std::string& format, const T&...args) {
	if(p >= Settings::getLogLevel()) {
		std::stringstream stream;
		printToStreamF(stream, format, args...);
		Impl::log(loc, p, stream.str());
	}
}


}} // namespace Aux::Log

#endif
