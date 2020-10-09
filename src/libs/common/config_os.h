#ifndef CONFIG_OS_H_
#define CONFIG_OS_H_


#define PESTPP_VERSION "5.0.3";

#if defined(_WIN32) || defined(_WIN64)
#define OS_WIN
#define DEF_DLAMCH DLAMCH
#define DEF_DLANBPRO_SPARCE DLANSVD
#define DEF_DLANSVD DLANSVD_SPARCE
#elif defined( __linux__)
#define OS_LINUX
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#define DEF_DLANSVD dlansvd_sparce_
#elif defined (__APPLE__)
#define OS_LINUX
#define DEF_DLAMCH dlamch_
#define DEF_DLANBPRO_SPARCE dlanbpro_sparce_
#define DEF_DLANSVD dlansvd_sparce_
#endif


#endif /* CONFIG_OS_H_ */
