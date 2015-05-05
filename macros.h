#ifndef MACROS_H
#define MACROS_H

#ifdef DEBUG
    #include <cassert>
    #define ASSERT(c) assert(c)
#else
    #define ASSERT(c)
#endif

#ifdef VERBOSE
    #define PRINT(c) std::cout << c << std::flush
#else
    #define PRINT(c)
#endif

#define NEXT(c) ((c+1)%3)
#define PREV(c) ((c+2)%3)

#endif // MACROS_H
