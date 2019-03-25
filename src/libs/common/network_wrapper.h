#ifndef NETWORK_H_
#define NETWORK_H_

#include "config_os.h"
#include "utilities.h"
//for windows
#ifdef OS_WIN
  #define FD_SETSIZE 2048
  #include <winsock2.h>
  #include <ws2tcpip.h>
  typedef int socklen_t;
#endif
//for linux
#ifdef OS_LINUX
  #include <arpa/inet.h>
  #include <unistd.h>
  #include <sys/types.h>
  #include <sys/select.h>
  #include<sys/wait.h>
  #include <errno.h>
  #include <signal.h>
  #include <netdb.h>
  #include <sys/socket.h>
#endif

//common for all systems
#include <iostream>
#include <vector>
#include <string>

void w_init();
int w_close(int sockfd);
void w_cleanup();
std::vector<std::string> w_getnameinfo_vec(int sockfd, int flags=0);
std::string w_getnameinfo_string(int sockfd, int flags = 0);
int w_getaddrinfo(const char *node, const char *service,
			  const struct addrinfo *hints, struct addrinfo **res);
int w_socket(int domain, int type, int protocol);
int w_connect(int sockfd, struct sockaddr *serv_addr, socklen_t addr_len);
int w_bind(int sockfd, struct sockaddr *my_addr, socklen_t addr_len);
int w_listen(int sockfd, int backlog);
int w_accept(int sockfd, struct sockaddr *addr, socklen_t *addr_len);
int w_send(int sockfd, int8_t *buf, int64_t len, int flags);
int w_sendall(int sockfd, int8_t *buf, int64_t *len);
int w_recv(int sockfd, int8_t *buf, int64_t len, int flags);
int w_recvall(int sockfd, int8_t *buf, int64_t *len);
int w_select(int numfds, fd_set *readfds, fd_set *writefds,
		   fd_set *exceptfds, struct timeval *timeout);
int w_memcpy_s(void *dest, size_t number_of_elements, const void *src, size_t count);
addrinfo* w_bind_first_avl(addrinfo *servinfo, int &sockfd);
addrinfo* w_connect_first_avl(addrinfo *servinfo, int &sockfd);
void w_print_servinfo(struct addrinfo *res, std::ostream &fout);
std::string w_get_addrinfo_string(struct addrinfo *p);
std::string w_get_error_msg();
void w_sleep(int millisec);
#endif /* NETWORK_H_ */

