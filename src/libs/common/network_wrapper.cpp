#include "network_wrapper.h"
#include "network_package.h"
#include "utilities.h"
#include "system_variables.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <thread>

#ifdef OS_WIN
#include <Windows.h>
#include <conio.h>
#endif


#ifdef OS_LINUX
#include <arpa/inet.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/select.h>
#include<sys/wait.h>
#include <errno.h>
#include <signal.h>
#endif

using namespace std;

string w_get_hostname()
{
	char* temp;
	string hostname("");
#ifdef OS_WIN
	temp = getenv("COMPUTERNAME");
	hostname = temp;
	temp = 0;

#else
	temp = getenv("HOSTNAME");
	if (temp != 0) {
		hostname = temp;
		temp = 0;
}
	else {
		temp = new char[512];
		if (gethostname(temp, 512) == 0) { // success = 0, failure = -1
			hostname = temp;
		}
		delete[]temp;
		temp = 0;
	}

#endif
	return hostname;

}

string w_init()
{
	stringstream ss;
#ifdef OS_WIN
	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2,0),  &wsaData)!=0)
	{
		ss << "WSAStartup failed " << w_get_error_msg() << endl;
	}
#endif
	return ss.str();
}

int w_close(int sockfd)
{
  int n;
	#ifdef OS_WIN
	shutdown(sockfd, SD_BOTH);
	if ((n = closesocket(sockfd)) != 0)
	{
		//cerr << "error closing socket: "  << w_get_error_msg() << endl;
	}
	return n;
	#endif
	#ifdef OS_LINUX
    n = shutdown(sockfd, 2);
	n = close(sockfd);
	#endif
	return n;
}

void w_cleanup()
{
   #ifdef OS_WIN
	WSACleanup();
   #endif
}


std::pair<int,std::string> w_getaddrinfo(const char *node, const char *service,
			  const struct addrinfo *hints, struct addrinfo **res)
{
	int status;
	stringstream ss;
	if ((status = getaddrinfo(node, service, hints, res)) !=0)
	{
		ss << "getaddrinfo error: " << gai_strerror(status) << endl;
	}
	return std::pair<int,string> (status,ss.str());
}

vector<string> w_getnameinfo_vec(int sockfd, int flags)
{
	int err;
	vector<string> name_info;
	char host[INET6_ADDRSTRLEN];
	char port[INET6_ADDRSTRLEN];
	struct sockaddr_storage addr;
	socklen_t addr_len = sizeof addr;
	err = getpeername(sockfd, (struct sockaddr*) &addr, &addr_len);
	err = getnameinfo((struct sockaddr*) &addr, addr_len, host, sizeof host, port, sizeof port, flags);
	name_info.push_back(host);
	name_info.push_back(port);
	return name_info;
}

string w_getnameinfo_string(int sockfd, int flags)
{
	vector<string> name_info_vec = w_getnameinfo_vec(sockfd, flags);
	stringstream ss;
	ss << name_info_vec[0] << ":" << name_info_vec[1];
	return ss.str();
}

int w_socket(int domain, int type, int protocol)
{
	int sockfd = socket(domain, type, protocol);
	if (sockfd < 0) {
		//cerr << "socket error: "  << w_get_error_msg() << endl;
	}
#ifdef OS_LINUX
//    int one = 1;
//    int retval = setsockopt(sockfd,SOL_SOCKET,SO_NOSIGPIPE,&one,sizeof(one));
//    if (retval != 0)
//    {
//        cout << "WARNING: unable to set SO_NOSIGPIPE on socket, continuing..." << endl;
//    }
    signal(SIGPIPE, SIG_IGN);

#endif
	return sockfd;
}

int w_connect(int sockfd, struct sockaddr *serv_addr, socklen_t addrlen)
{
	int n=0;
	if ((n=connect(sockfd, serv_addr, addrlen)) == -1 )
	{
		//cerr << "connect error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_bind(int sockfd, struct sockaddr *my_addr, socklen_t addrlen)
{
	int n=0;
	if ((n=::bind(sockfd, my_addr, addrlen)) == -1 )
	{

		//cerr << "bind error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_accept(int sockfd, struct sockaddr *addr, socklen_t *addrlen)
{
	int n=0;
	if ((n=accept(sockfd, addr, addrlen)) == -1)
	{
		//cerr << "bind error: " << w_get_error_msg() << endl;
	}
	return n;
}


int w_listen(int sockfd, int backlog)
{
	int n;
	if ((n = listen(sockfd, backlog)) == -1)
	{
		//cerr << "listen error: "  << w_get_error_msg() << endl;
	}
	return n;
}

int w_recv(int sockfd, int8_t *buf, int64_t len, int flags)
{
	int n;
	n = recv(sockfd, (char*)buf, len, flags);
	if (n < 0){
		//cerr << "recv error: "  << w_get_error_msg() << endl;
	}
	return n;
}
int w_send(int sockfd, int8_t *buf, int64_t len, int flags)
{
	int n;
	n = send(sockfd, (char*)buf, len, flags);
	if (n < 0){
		//cerr << "send error: "  << w_get_error_msg() << endl;
	}
	return n;
}

int w_sendall(int sockfd, int8_t *buf, int64_t *len)
{
	unsigned long total = 0; // how many bytes we've sent
	unsigned long bytesleft = *len; // how many we have left to send
	int n;
	while(total < *len) {
		n = send(sockfd, (char*)buf + total, bytesleft, 0);
		if (n == -1) { break; }  //error
		if (n == 0) { break; } //connection closed
		total += n;
		bytesleft -= n;
	}
	*len = total; // return number actually sent here
		if (n < 0){
		//cerr << "w_sendall error: " << n << endl;
	}
	if (n > 0) {n = 1;}
	return n; // return -1 on failure, 0 closed connection or 1 on success
}


int w_recvall(int sockfd, int8_t *buf, int64_t *len)
{
	unsigned long total = 0; // how many bytes we've received
	unsigned long bytesleft = *len; // how many we have left to receive
	int n;
	while(total < *len) {
		n = recv(sockfd, (char*)buf + total, bytesleft, 0);
		if (n == -1) { break; }  //error
		if (n == 0) { break; } //connection closed
		total += n;
		bytesleft -= n;
	}
	*len = total; // return number actually received here
	if (n < 0){n = -1;}
	else if (n > 0) {n = 1;}
	return n; // return -1 on failure, 0 closed connection or 1 on success
}

addrinfo* w_bind_first_avl(addrinfo *servinfo, int &sockfd)
{
	// loop through all the results and bind to the first we can
	struct addrinfo *p;
	char yes = '1';
	for(p = servinfo; p != nullptr; p = p->ai_next)
	{
		if ((sockfd = w_socket(p->ai_family, p->ai_socktype,
		p->ai_protocol)) == -1)
		{
			continue;
		}
		if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes,
		sizeof(int)) == -1)
		{
			return nullptr;
		}
		if (w_bind(sockfd, p->ai_addr, p->ai_addrlen) == -1)
		{
			w_close(sockfd);
			continue;
		}
	break;
	}
	if (p == nullptr) {
		//cerr << "server: failed to bind" << endl;
	}
	return p;
}

addrinfo* w_connect_first_avl(addrinfo *servinfo, int &sockfd)
{
	// loop through all the results and connect to the first we can
	struct addrinfo *p;
	for(p = servinfo; p != nullptr; p = p->ai_next)
	{
		if ((sockfd = w_socket(p->ai_family, p->ai_socktype,
		p->ai_protocol)) == -1)
		{
			continue;
		}
		if (w_connect(sockfd, p->ai_addr, p->ai_addrlen) == -1)
		{
			w_close(sockfd);
			continue;
		}
	break;
	}
	if (p == nullptr) {
		// connection failed
	}
	return p;
}



void w_print_servinfo(addrinfo *res, ostream &fout)
{
	struct addrinfo *p;
	fout << "IP addresses:" << endl;
	for(p = res;p != NULL; p = p->ai_next)
	{
		string socket_string = w_get_addrinfo_string(p);
	   fout << "  " << socket_string << endl;
	}
}

string w_get_addrinfo_string(struct addrinfo *p)
{
	stringstream sstr;
	void *addr;
	
	string ipver;
	unsigned short port = 0;
	// get the pointer to the address itself,
	// different fields in IPv4 and IPv6:
	if (p->ai_family == AF_INET) { // IPv4
		struct sockaddr_in *ipv4 = (struct sockaddr_in *)p->ai_addr;
		addr = &(ipv4->sin_addr);
		port = ntohs(ipv4->sin_port);
		ipver = "IPv4";
		char ipstr[INET_ADDRSTRLEN];
		inet_ntop(p->ai_family, addr, ipstr, sizeof ipstr);
		sstr << ipstr <<":" << port<< " (" << ipver << ")";
	}
	else { // IPv6
		struct sockaddr_in6 *ipv6 = (struct sockaddr_in6 *)p->ai_addr;
		addr = &(ipv6->sin6_addr);
		port = ntohs(ipv6->sin6_port);
		ipver = "IPv6";
		char ipstr[INET6_ADDRSTRLEN];
		inet_ntop(p->ai_family, addr, ipstr, sizeof ipstr);
		sstr << "[" << ipstr <<"]:" << port << " (" << ipver << ")";
	}
	// convert the IP to a string and print it:
	return sstr.str();
}


int w_select(int numfds, fd_set *readfds, fd_set *writefds,
		   fd_set *exceptfds, struct timeval *timeout)
{
	int n;
	if ((n=select(numfds, readfds, writefds, exceptfds, timeout)) == -1)
	{
		cerr << "select error: " << w_get_error_msg() << endl;
	}
	return n;
}

int w_memcpy_s(void *dest, size_t numberOfElements, const void *src, size_t count)
{
	int err = 0;
		#ifdef OS_WIN
	err = memcpy_s(dest, numberOfElements, src, count);
		#endif
		#ifdef OS_LINUX
	memcpy(dest, src, count);
		#endif
	if (err) {
		cerr << "Error executing memcpy" << endl;
	}
	return err;
}

string w_get_error_msg()
{
	stringstream err_msg;
	#ifdef OS_WIN
	int err_no = WSAGetLastError();
	LPVOID lpMsgBuf;
	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER |
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		err_no,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
		(LPTSTR) &lpMsgBuf,
		0,
		NULL);
	// Process any inserts in lpMsgBuf.
	// ..
	err_msg << "(tcp/ip error id = " << errno << ") " << (LPCTSTR)lpMsgBuf ;
	// Free the buffer.
	LocalFree( lpMsgBuf );
	#endif
	#ifdef OS_LINUX
	err_msg << "(tcp/ip error id = " << errno << ") " << strerror(errno);
	#endif
	return err_msg.str();
}

void w_sleep(int millisec)
{
   #ifdef OS_WIN
	Sleep(millisec);
   #endif
   #ifdef OS_LINUX
	sleep(millisec / 1000);
   #endif
}

