#ifndef CONFIG_OS_H_
#define CONFIG_OS_H_
#include <string>

#define PESTPP_VERSION "5.2.24";


#if defined(_WIN32) || defined(_WIN64)
#define OS_WIN
#define DEF_DLAMCH DLAMCH
#define DEF_DLANBPRO_SPARCE DLANSVD
#define DEF_DLANSVD DLANSVD_SPARCE
#define OS_SEP '\\'
#elif defined (__linux__)
#define OS_LINUX
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#define DEF_DLANSVD dlansvd_sparce_
#define OS_SEP '/'
#elif defined (__APPLE__)
#define OS_LINUX
#define OS_MAC
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#define DEF_DLANSVD dlansvd_sparce_
#define OS_SEP '/'
#endif


#endif /* CONFIG_OS_H_ */
