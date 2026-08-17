#ifndef CONFIG_OS_H_
#define CONFIG_OS_H_
#include <string>

// The one place a version is defined. cmake parses this line to derive the package version,
// so editing it here is the whole version bump.
//
// A release candidate carries a suffix directly on the string - "5.2.29rc1" - and a final
// release has none. cmake splits the numeric major.minor.patch off for project() and CPack's
// numeric fields, and uses the FULL string for the package name, so an rc build produces
// pestpp-5.2.29rc1-mac.tar.gz and cannot be mistaken for the release.
//
// No trailing semicolon. There used to be one, which made this macro expand to `"5.2.28";` -
// so it happened to work as `string v = PESTPP_VERSION;` and would not compile anywhere else,
// string(PESTPP_VERSION) included.
#define PESTPP_VERSION "6.0.0a"

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
