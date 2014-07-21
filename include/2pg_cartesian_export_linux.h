
#ifndef _2PG_CARTESIAN_EXPORT_H
#define _2PG_CARTESIAN_EXPORT_H

#ifdef _2PG_CARTESIAN_STATIC_DEFINE
#  define _2PG_CARTESIAN_EXPORT
#  define _2PG_CARTESIAN_NO_EXPORT
#else
#  ifndef _2PG_CARTESIAN_EXPORT
#    ifdef _2pg_cartesian_EXPORTS
        /* We are building this library */
#      define _2PG_CARTESIAN_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define _2PG_CARTESIAN_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef _2PG_CARTESIAN_NO_EXPORT
#    define _2PG_CARTESIAN_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef _2PG_CARTESIAN_DEPRECATED
#  define _2PG_CARTESIAN_DEPRECATED __attribute__ ((__deprecated__))
#  define _2PG_CARTESIAN_DEPRECATED_EXPORT _2PG_CARTESIAN_EXPORT __attribute__ ((__deprecated__))
#  define _2PG_CARTESIAN_DEPRECATED_NO_EXPORT _2PG_CARTESIAN_NO_EXPORT __attribute__ ((__deprecated__))
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define _2PG_CARTESIAN_NO_DEPRECATED
#endif

#endif
