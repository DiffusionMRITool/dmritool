// smart_assert.h
//
//////////////////////////////////////////////////////////////////////

/** 
 *  The code is modified from smart_assert by John Torjo 
 *  http://www.torjo.com 
 *  http://www.drdobbs.com/cpp/enhancing-assertions/184403745 
 *
 *  \author  John Torjo, Jian Cheng
 * */

#if !defined(SMART_ASSERT_H)
#define SMART_ASSERT_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if _MSC_VER > 1000

// note:
// moving this after pragma push will render it useless (VC6)
//
// identifier truncated to 255 chars in debug information
#pragma warning ( disable : 4786)

#pragma warning ( push )
// *this used in base-member initialization; it's ok
#pragma warning ( disable : 4355)
#endif

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <map>

#include <fstream>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <stdexcept>

#include "utlCoreMacro.h"


#ifdef MATLAB_MEX_FILE
#include <mex.h>
#define SA_Abort(expout) mexErrMsgTxt(expout)
#else
#define SA_Abort(expout) do { std::cerr << expout <<"\n" << std::flush; abort(); } while(0)
#endif

#if defined(__BORLANDC__)
  #define __SMART_ASSERT_LOCATION__ __FUNC__
#elif defined(_WIN32) && !defined(__MINGW32__) && !defined(__CYGWIN__) && !defined(CSWIG)
  #define __SMART_ASSERT_LOCATION__ __FUNCSIG__
#elif defined(__GNUC__)
  #define __SMART_ASSERT_LOCATION__ __PRETTY_FUNCTION__
#else
  #define __SMART_ASSERT_LOCATION__ __FUNCTION__
#endif


namespace smart_assert {
enum {

    // default behavior - just loggs this assert, 
    lvl_log = 0,

    // default behavior - just loggs this assert
    // (a message is shown to the user to the console)
    lvl_warn = 100,

    // default behavior - asks the user what to do:
    // Ignore/ Retry/ etc.
    lvl_debug = 200,

    // default behavior - throws a smart_assert_error
    lvl_error = 300,

    // default behavior - dumps all assert context to console, 
    // and aborts
    lvl_fatal = 1000
};

enum{
  lvl_condition_assert = 0, 
  lvl_condition_exception = 1
};


/* 
    contains details about a failed assertion
*/
class assert_context {
    typedef std::string string;
public:
    assert_context() : level_( lvl_log), condition_(lvl_condition_assert) {
    }

    // where the assertion failed: file & line
    void set_file_line( const char * file, int line) {
        file_ = file;
        line_ = line;
    }
    void set_file_line_func_cond( const char * file, int line, const char * func, int cond) {
        file_ = file;
        line_ = line;
        func_ = func;
        condition_ = cond;
    }
    const string & get_context_file() const { return file_; }
    const string & get_context_func() const { return func_; }
    int get_context_line() const { return line_; }

    // get/ set expression
    void set_expr( const string & str) { expr_ = str; }
    const string & get_expr() const { return expr_; }

    typedef std::pair< string, string> val_and_str;
    typedef std::vector< val_and_str> vals_array;
    // return values array as a vector of pairs:
    // [Value, corresponding string]
    const vals_array & get_vals_array() const { return vals_; }
    // adds one value and its corresponding string
    void add_val( const string & val, const string & str) {
        vals_.push_back( val_and_str( val, str) );
    }

    // get/set level of assertion
    void set_level( int nLevel) { level_ = nLevel; }
    int get_level() const { return level_; }
    
    void set_condition( int cond) { condition_ = cond; }
    int get_condition() const { return condition_; }

    // get/set (user-friendly) message 
    void set_level_msg( const char * strMsg)  { 
        if ( strMsg)
            msg_ = strMsg; 
        else
            msg_.erase();
    }
    const string & get_level_msg() const { return msg_; }

private:
    // where the assertion occured
    string file_;
    string func_;
    int line_;

    // expression and values
    string expr_;
    vals_array vals_;

    // level and message
    int level_;
    int condition_;
    string msg_;
};



    typedef void (*assert_func)( const assert_context & context);

    // helpers
    inline std::string get_typeof_level( int nLevel);
    inline void dump_context_summary( const assert_context & context, std::ostream & out);
    inline void dump_context_detail( const assert_context & context, std::ostream & out);

    // defaults
    inline void default_log_handler( const assert_context & context);
    inline void default_warn_handler( const assert_context & context);
    inline void default_debug_handler( const assert_context & context);
    inline void default_error_handler( const assert_context & context);
    inline void default_fatal_handler( const assert_context & context);
    inline void default_logger( const assert_context & context);


// namespace Private {
    inline void init_assert();
    inline void set_default_log_stream( std::ostream & out);
    inline void set_default_log_name( const char * str);

    // allows finding if a value is of type 'const char *'
    // and is null; if so, we cannot print it to an ostream
    // directly!!!
    template< class T>
    struct is_null_finder {
        bool is( const T &) const {
            return false;
        }
    };

    template<>
    struct is_null_finder< char*> {
        bool is( char * const & val) {
            return val == 0;
        }
    };

    template<>
    struct is_null_finder< const char*> {
        bool is( const char * const & val) {
            return val == 0;
        }
    };


// } // namespace Private


struct Assert {
    typedef smart_assert::assert_func assert_func;

    // helpers, in order to be able to compile the code
    Assert & SMART_ASSERT_A;
    Assert & SMART_ASSERT_B;

    Assert( const char * expr) 
        : SMART_ASSERT_A( *this), 
          SMART_ASSERT_B( *this),
          needs_handling_( true) {
        context_.set_expr( expr);

        if ( ( logger() == 0) || handlers().size() < 5) {
            // used before main!
            init_assert();
        }
    }

    Assert( const Assert & other)
        : SMART_ASSERT_A( *this), 
          SMART_ASSERT_B( *this),
          context_( other.context_),
          needs_handling_( true) {
        other.needs_handling_ = false;
    }

    ~Assert() {
        if ( needs_handling_) 
            handle_assert();
    }

    template< class type>
    Assert & print_current_val( const type & val, const char * msg) {
        std::ostringstream out;

        is_null_finder< type> f;
        bool bIsNull = f.is( val);
        if ( !bIsNull)
            out << val;
        else
            // null string
            out << "null";
        context_.add_val( out.str(), msg);
        return *this;
    }

    Assert & print_context( const char * file, int line) {
        context_.set_file_line( file, line);
        return *this;
    }
    Assert & print_context( const char * file, int line, const char* func, int cond) {
        context_.set_file_line_func_cond( file, line, func, cond);
        return *this;
    }

    Assert & msg( const char * strMsg) {
        context_.set_level_msg( strMsg);
        return *this;
    }

    Assert & level( int nLevel, const char * strMsg = 0) {
        context_.set_level( nLevel);
        context_.set_level_msg( strMsg);
        return *this;
    }

    Assert & log( const char * strMsg = 0) {
        return level( lvl_log, strMsg);
    }

    Assert & warn( const char * strMsg = 0) {
        return level( lvl_warn, strMsg);
    }

    Assert & debug( const char * strMsg = 0) {
        return level( lvl_debug, strMsg);
    }

    Assert & error( const char * strMsg = 0) {
        return level( lvl_error, strMsg);
    }

    Assert & fatal( const char * strMsg = 0) {
        return  level( lvl_fatal, strMsg);
    }

    // in this case, we set the default logger, and make it
    // write everything to this file
    static void set_log( const char * strFileName) {
        set_default_log_name( strFileName);
        logger() = &smart_assert::default_logger;
    }

    // in this case, we set the default logger, and make it
    // write everything to this log
    static void set_log( std::ostream & out) {
        set_default_log_stream( out);
        logger() = &smart_assert::default_logger;
    }

    static void set_log( assert_func log) {
        logger() = log;
    }

    static void set_handler( int nLevel, assert_func handler) {
        handlers()[ nLevel] = handler;
    }

private:
    // handles the current assertion.
    void handle_assert() {
        logger()( context_);
        get_handler( context_.get_level() )( context_);
    }

    /*
        IMPORTANT NOTE:
        The only reason logger & handlers are functions, are
        because you might use SMART_ASSERT before main().

        In this case, since they're statics, they might not
        be initialized. However, making them functions 
        will make it work.
    */

    // the log
    static assert_func & logger() {
        static assert_func inst;
        return inst;
    }

    // the handler
    typedef std::map< int, assert_func> handlers_collection; 
    static handlers_collection & handlers() {
        static handlers_collection inst;
        return inst;
    }

    static assert_func get_handler( int nLevel) {
        handlers_collection::const_iterator found = handlers().find( nLevel);
        if ( found != handlers().end() )
            return found->second;
        else
            // we always assume the debug handler has been set
            return handlers().find( lvl_debug)->second;
    }

private:
    assert_context context_;
    mutable bool needs_handling_;

};

    inline Assert make_assert( const char * expr) {
        return Assert( expr);
    }
} // namespace smart_assert



////////////////////////////////////////////////////////
// macro trickery

// note: NEVER define SMART_ASSERT_DEBUG directly
// (it will be overridden);
//
// #define SMART_ASSERT_DEBUG_MODE instead

#ifdef SMART_ASSERT_DEBUG_MODE 
    #if SMART_ASSERT_DEBUG_MODE == 1
    #define SMART_ASSERT_DEBUG 
    #else
    #undef SMART_ASSERT_DEBUG
    #endif

#else

// defaults
    #ifndef NDEBUG
    #define SMART_ASSERT_DEBUG 
    #else
    #undef SMART_ASSERT_DEBUG
    #endif
#endif

#ifdef SMART_ASSERT_DEBUG
// "debug" mode
#define SMART_ASSERT( expr) \
    if ( (expr) ) ; \
    else ::smart_assert::make_assert( #expr).print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_assert).SMART_ASSERT_A \
    /**/

#else
// "release" mode
#define SMART_ASSERT( expr) \
    if ( true ) ; \
    else ::smart_assert::make_assert("").SMART_ASSERT_A \
    /**/

#endif // ifdef SMART_ASSERT_DEBUG


#define SMART_VERIFY( expr) \
    if ( (expr) ) ; \
    else ::smart_assert::make_assert( #expr).error().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_assert).SMART_ASSERT_A \
    /**/

#define SMART_EXCEPTION( expr) \
    if ( !(expr) ) ; \
    else ::smart_assert::make_assert( #expr).error().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A \

#define SMART_PRINT \
    if ( false ) ; \
    else ::smart_assert::make_assert("").log().print_context( __FILE__, __LINE__,__SMART_ASSERT_LOCATION__, ::smart_assert::lvl_condition_exception).SMART_ASSERT_A \

#define SMART_ASSERT_A(x) SMART_ASSERT_OP(x, B)
#define SMART_ASSERT_B(x) SMART_ASSERT_OP(x, A)

#define SMART_ASSERT_OP(x, next) \
    SMART_ASSERT_A.print_current_val((x), #x).SMART_ASSERT_##next \
    /**/


#if _MSC_VER > 1000
#pragma warning ( pop )
#endif

inline void break_into_debugger() {
// MSVC, BCB, 
#if (defined _MSC_VER) || (defined __BORLANDC__)
        __asm { int 3 };
#elif defined(__GNUC__)
    // GCC
    __asm ("int $0x3");
#else
    #  error Please supply instruction to break into code
#endif
}


namespace {
    // in case we're logging using the default logger...
    struct stream_holder {
        stream_holder() : out_( 0), owns_( false) {}
        ~stream_holder() {
            if ( owns_)
                delete out_;
            out_ = 0;
        }
        std::ostream * out_;
        bool owns_;
    };
    // information about the stream we write to, in case 
    // we're using the default logger
    stream_holder default_logger_info;

    // intitializes the SMART_ASSERT library
    struct assert_initializer {
        assert_initializer() {
          smart_assert::init_assert();
        }
    } init;
} // anonymous namespace



namespace smart_assert {

    // returns a message corresponding to the type of level
    inline std::string get_typeof_level( int nLevel) {
        switch ( nLevel) {
        case lvl_log: return  __UTL_LOG_STRING;
        case lvl_warn: return __UTL_WARNING_STRING;
        case lvl_debug: return __UTL_DEBUG_STRING;
        case lvl_error: return __UTL_ERROR_STRING;
        case lvl_fatal: return __UTL_FATAL_STRING;
        default: {
            std::ostringstream out;
            out << "(level=" << nLevel << ")";
            return out.str();
                 }
        };
    }

    // helpers, for dumping the assertion context
    inline void dump_context_summary( const assert_context & context, std::ostream & out) {
        out << "\n" << get_typeof_level( context.get_level() ) 
            << " in "<<__UTL_BOLD("File")<<": " << context.get_context_file() << ", "<<__UTL_BOLD("Line")<<": " << context.get_context_line() << ", "<<__UTL_BOLD("Function")<<": " << context.get_context_func() << '\n';
        if ( !context.get_level_msg().empty())
            // we have a user-friendly message
            out << context.get_level_msg();
        else
          if (context.get_expr()!="\"\"" && context.get_expr()!="")
            out << __UTL_BOLD("Expression")<<" : '" << __UTL_EXPSTR(context.get_expr()) <<"' " << (context.get_condition()==lvl_condition_assert?"failed":"satisfied");
        out << std::endl;
    }

    inline void dump_context_detail( const assert_context & context, std::ostream & out) {
        out << "\n" << get_typeof_level( context.get_level() ) 
            << " in "<<__UTL_BOLD("File")<<": " << context.get_context_file() << ", "<<__UTL_BOLD("Line")<<": " << context.get_context_line() << ", "<<__UTL_BOLD("Function")<<": " << context.get_context_func() << '\n';
        if ( !context.get_level_msg().empty())
            out << __UTL_BOLD("msg")<<": '" << context.get_level_msg() << "'\n";
        if (context.get_expr()!="\"\"" && context.get_expr()!="")
          out << __UTL_BOLD("Expression")<<" : '" << __UTL_EXPSTR(context.get_expr()) <<"' " << (context.get_condition()==lvl_condition_assert?"failed":"satisfied") <<"\n";
        
        typedef assert_context::vals_array vals_array;
        const vals_array & aVals = context.get_vals_array();
        if ( !aVals.empty() ) {
            bool bFirstTime = true;
            vals_array::const_iterator first = aVals.begin(), last = aVals.end();
            while ( first != last) {
                if ( bFirstTime) {
                    out << "Values: ";
                    bFirstTime = false;
                }
                else {
                    out << "        ";
                }
                out << first->second << "='" << first->first << "'\n";
                ++first;
            }
        }
        out << std::endl;
    }
    
    inline void dump_context_log_detail( const assert_context & context, std::ostream & out) {
        out << "\n" << get_typeof_level( context.get_level() ) 
            << " in "<<__UTL_BOLD("File")<<": " << context.get_context_file() << ", "<<__UTL_BOLD("Line")<<": " << context.get_context_line() << ", "<<__UTL_BOLD("Function")<<": " << context.get_context_func() << '\n';
        if ( !context.get_level_msg().empty())
            out << __UTL_BOLD("msg")<<": '" << context.get_level_msg() << "'\n";
        if (context.get_expr()!="\"\"" && context.get_expr()!="")
          out << __UTL_BOLD("Expression")<<" : '" << __UTL_EXPSTR(context.get_expr()) <<"' " << (context.get_condition()==lvl_condition_assert?"failed":"satisfied") <<"\n";
        
        typedef assert_context::vals_array vals_array;
        const vals_array & aVals = context.get_vals_array();
        if ( !aVals.empty() ) {
          if (aVals.size()==1)
            out << "(" << aVals[0].second << ") = " << "(" << aVals[0].first <<  ")"; 
          else
            {
            out << "(";
            for ( int i = 0; i < aVals.size()-1; i += 1 ) 
              out << aVals[i].second << ", ";
            out << aVals.back().second << ") = (";
            for ( int i = 0; i < aVals.size()-1; i += 1 ) 
              out << aVals[i].first << ", ";
            out << aVals.back().first << ")";
            }
          out << std::endl;
        }
        out << std::endl;
    }

    ///////////////////////////////////////////////////////
    // logger

    inline void default_logger( const assert_context & context) {
        if ( default_logger_info.out_ == 0)
            return;
        dump_context_log_detail( context, *( default_logger_info.out_) );
      // if ( default_logger_info.out_ == 0)
      //   dump_context_log_detail( context, std::cout );
      // else
      //   dump_context_log_detail( context, *( default_logger_info.out_));
    }

    ///////////////////////////////////////////////////////
    // handlers
    
    inline void default_log_handler( const assert_context & context)
      {
      if ( default_logger_info.out_ == 0)
        dump_context_log_detail( context, std::cout );
      else
        dump_context_log_detail( context, *( default_logger_info.out_));
      }

    // warn : just dump summary to console
    inline void default_warn_handler( const assert_context & context) {
        dump_context_summary( context, std::cout);
    }


    // debug: ask user what to do
   inline void default_debug_handler( const assert_context & context) {
        static bool ignore_all = false;
        if ( ignore_all)
            // ignore All asserts
            return;
        typedef std::pair< std::string, int> file_and_line;
        static std::set< file_and_line> ignorer;
        if ( ignorer.find( file_and_line( context.get_context_file(), context.get_context_line())) != ignorer.end() )
            // this is Ignored Forever
            return;

        dump_context_summary( context, std::cerr );
        std::cerr << "\nPress (I)gnore/ Igore (F)orever/ Ignore (A)ll/ (D)ebug/ A(b)ort: ";
        std::cerr.flush();
        char ch = 0;

        bool bContinue = true;
        while ( bContinue && std::cin.get( ch)) {
            bContinue = false;
            switch ( ch) {
            case 'i': case 'I':
                // ignore
                break;

            case 'f': case 'F':
                // ignore forever
                ignorer.insert( file_and_line( context.get_context_file(), context.get_context_line()));
                break;

            case 'a': case 'A':
                // ignore all
                ignore_all = true;
                break;

            case 'd': case 'D':
                // break
                break_into_debugger();
                break;

            case 'b': case 'B':
                SA_Abort("");
                break;

            default:
                bContinue = true;
                break;
            }
        }
    }


    // error : throw a runtime exception
    inline void default_error_handler( const assert_context & context) {
        std::ostringstream out;
        dump_context_detail( context, out);
        throw std::runtime_error( out.str());
    }


    // fatal : dump error and abort
    inline void default_fatal_handler( const assert_context & context) {
        dump_context_detail( context, std::cerr);
        SA_Abort("");
    }





    inline void init_assert() {
        Assert::set_log( &::smart_assert::default_logger);
        Assert::set_handler( lvl_log, &::smart_assert::default_log_handler);
        Assert::set_handler( lvl_warn, &::smart_assert::default_warn_handler);
        Assert::set_handler( lvl_debug, &::smart_assert::default_debug_handler);
        Assert::set_handler( lvl_error, &::smart_assert::default_error_handler);
        Assert::set_handler( lvl_fatal, &::smart_assert::default_fatal_handler);
    }

    // sets the default logger to write to this stream
    inline void set_default_log_stream( std::ostream & out) {
        default_logger_info.out_ = &out;
        default_logger_info.owns_ = false;
    }

    // sets the default logger to write to this file
    inline void set_default_log_name( const char * str) {
        default_logger_info.owns_ = false;
        default_logger_info.out_ = new std::ofstream( str);
        default_logger_info.owns_ = true;
    }

} // namespace smart_assert


#endif 
