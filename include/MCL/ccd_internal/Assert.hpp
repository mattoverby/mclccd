// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_ASSERT_HPP
#define MCL_ASSERT_HPP 1

#include <cstdlib>
#include <string>
#include <stdexcept>

namespace mcl
{

static inline void mclAssertHandler(bool cond, const std::string& file, const int& line)
{
    if (!cond)
    {
        std::string err_msg = "Assertion failed in "+file+" line "+std::to_string(line);
        throw std::runtime_error(err_msg.c_str());
    }
}

static inline void mclAssertHandlerMsg(bool cond, const std::string& file, const int& line, const std::string &msg)
{
    if (!cond)
    {
        std::string err_msg = "Assertion failed in "+file+" line "+std::to_string(line)+": "+msg;
        throw std::runtime_error(err_msg.c_str());
    }
}

} // ns mcl

// Neat trick that allows macros with multiple arguments:
// https://stackoverflow.com/questions/3046889/optional-parameters-with-c-macros
// Also consider using #condition

#define mclAssert_withmsg(cond, msg) mcl::mclAssertHandlerMsg(cond, std::string(__FILE__), __LINE__, msg)
#define mclAssert_nomsg(cond) mcl::mclAssertHandler(cond, std::string(__FILE__), __LINE__)
#define mclAssert_stripargs(xx,cond,msg,FUNC, ...) FUNC

#define mclAssert(...) \
    mclAssert_stripargs(,##__VA_ARGS__,\
    mclAssert_withmsg(__VA_ARGS__),\
    mclAssert_nomsg(__VA_ARGS__),\
    )

#endif