#include <memory>
#include <sstream>
#include <cstring>
#include "network_package.h"
#include "network_wrapper.h"
#include <cassert>

using namespace std;

//Static Memeber Initialization
int64_t NetPackage::last_group_id = 0;
int8_t NetPackage::security_code[5] = { 1, 5, 25, 50, 100 };

//Static Methods
int NetPackage::get_new_group_id()
{
	return ++last_group_id;
}

bool NetPackage::allowable_ascii_char(int8_t value)
{
	// This is for security resasons.  Only \0
	//  and printable characters are allowed
	bool ret_val = false;
	if (value == 0 || (value >= 32 && value <= 176))
		ret_val = true;
	return ret_val;
}

bool NetPackage::check_string(const int8_t *data_src, size_t _size)
{
	bool safe_data = true;
	for (size_t i = 0; i < _size; ++i)
	{
		if (!allowable_ascii_char(data_src[i]))
		{
			safe_data = false;
			return safe_data;
		}
	}
	return safe_data;
}

bool NetPackage::check_string(const vector<int8_t> &data_src, size_t index1, size_t _size)
{
	size_t index2 = min(index1 + _size, data_src.size());
	bool safe_data = check_string(&(data_src[index1]), index2 - index1);
	return safe_data;
}

string NetPackage::extract_string(const int8_t *data_src, size_t _size)
{
	vector<char> buf;
	// This is done to remove possible system dependicies on whether char/uchar
	// is use to represent a standard char
	for (size_t i = 0; i < _size; ++i)
	{
		buf.push_back(data_src[i]);
	}
	string ret_val(buf.begin(), buf.end());
	return ret_val;
}

string NetPackage::extract_string(const vector<int8_t> &data_src, size_t index1, size_t _size)
{
	size_t index2 = min(index1 + _size, data_src.size());
	string  ret_val = extract_string(&(data_src[index1]), index2 - index1);
	return ret_val;
}

//template<class InputIterator>
//string NetPackage::extract_string(InputIterator first, InputIterator last)
//{
//	vector<char> buf;
//	while (first != last) {
//		buf.push_back(*first);
//		++first;
//	}
//	string ret_val(buf.begin(), buf.end());
//	return ret_val;
//}

template<class InputIterator>
vector<int8_t> NetPackage::pack_string(InputIterator first, InputIterator last)
{
	vector<int8_t> buf;
	while (first != last) {
		buf.push_back(*first);
		++first;
	}
	return buf;
}

//Non static methods
NetPackage::NetPackage(PackType _type, int _group, int _run_id, const string &desc_str)
	: type(_type), group(_group), run_id(_run_id)
{
	memset(desc, '\0', DESC_LEN);
	int i = 0;
	int max_len = min(size_t(DESC_LEN - 1), desc_str.size());
	// This is done to remove possible system dependicies on whether char/uchar
	// is use to represent a standard char
	int i_desc = 0;
	for (int i = 0; i < max_len; ++i)
	{
		if (!allowable_ascii_char(desc_str[i]))
		{
		}
		else
		{
			desc[i_desc++] = desc_str[i];
		}
	}
	//strncpy(desc, desc_str.c_str(), DESC_LEN-1);
	data_len = 1;
}
void NetPackage::reset(PackType _type, int _group, int _run_id, const string &_desc)
{
	type = _type;
	group = _group;
	run_id = _run_id;
	memset(desc, '\0', DESC_LEN);
	int max_len = min(size_t(DESC_LEN - 1), _desc.size());
	// This is done to remove possible system dependicies on whether char/uchar
	// is use to represent a standard char
	int i_desc = 0;
	for (int i = 0; i < max_len; ++i)
	{
		if (!allowable_ascii_char(_desc[i]))
		{
		}
		else
		{
			desc[i_desc++] = _desc[i];
		}
	}
	//strncpy(desc, _desc.c_str(), DESC_LEN-1);
	data.clear();
}

int NetPackage::send(int sockfd, const void *data, int64_t data_len_l)
{
	int n;

	// first send security code
	int64_t security_code_size = sizeof(security_code);
	n = w_sendall(sockfd, &security_code[0], &security_code_size);
	if (n<1) {
		cerr << "NetPackage::send error: could not send security code" << endl;
		return n;
	}
	int64_t buf_sz = 0;
	//calculate the size of buffer
	buf_sz += sizeof(buf_sz);
	buf_sz += sizeof(type);
	buf_sz += sizeof(group);
	buf_sz += sizeof(run_id);
	buf_sz += sizeof(desc);
	buf_sz += data_len_l;
	//pack information into buffer
	//unique_ptr<char[]> buf(new char[buf_sz]);
	vector<int8_t> buf;
	buf.resize(buf_sz, '\0');
	size_t i_start = 0;
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &buf_sz, sizeof(buf_sz));
	i_start += sizeof(buf_sz);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &type, sizeof(type));
	i_start += sizeof(type);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &group, sizeof(group));
	i_start += sizeof(group);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, &run_id, sizeof(run_id));
	i_start += sizeof(run_id);
	w_memcpy_s(&buf[i_start], buf_sz-i_start, desc, sizeof(desc));
	i_start += sizeof(desc);
	if (data_len_l > 0) {
		w_memcpy_s(&buf[i_start], buf_sz-i_start, data, data_len_l);
		i_start += data_len_l;
	}
	n = w_sendall(sockfd, buf.data(), &buf_sz);
	if (i_start != buf_sz) {
		cerr << "NetPackage::send error: could only send" << i_start
			<< " out of " << buf_sz << "bytes" << endl;
		n = -2;
	}
	return n;  // return -2 on corrupt send, -1 on failure, 0 closed connection or 1 on success
}

int  NetPackage::recv(int sockfd)
{
	long n;
	int64_t header_sz = 0;
	int64_t buf_sz = 0;
	size_t i_start = 0;
	int8_t rcv_security_code[5] = { 0, 0, 0, 0, 0 };
	int64_t rcv_security_code_size = sizeof(rcv_security_code);
	bool corrupt_desc = false;

	try{
		//get header (ie size, seq_id, id and name)
		header_sz = sizeof(buf_sz) + sizeof(type) + sizeof(group) + sizeof(run_id) + sizeof(desc);
		vector<int8_t> header_buf;
		header_buf.resize(header_sz, '\0');
		n = w_recvall(sockfd, &rcv_security_code[0], &rcv_security_code_size);
		int security_cmp = memcmp(security_code, rcv_security_code, sizeof(security_code));
		if (security_cmp != 0)
		{
			// corrupt message; message did not originate from a PEST++ application
			n = -2;
			cerr << "NetPackage::recv error - message did not originate from a PEST++ application" << endl;
			return n;
		}

		n = w_recvall(sockfd, &header_buf[0], &header_sz);
		if (n > 0 && header_sz != header_buf.size()) {
			// corrupt message; message not the correct length
			n = -2;
			cerr << "NetPackage::recv error reading header: expected" << header_buf.size()
				<< " bytes, but received " << header_sz << "bytes" << endl;
		}
		else if (n > 0) {
			i_start = 0;
			w_memcpy_s(&buf_sz, sizeof(buf_sz), &header_buf[i_start], sizeof(buf_sz));
			i_start += sizeof(buf_sz);
			w_memcpy_s(&type, sizeof(type), &header_buf[i_start], sizeof(type));
			i_start += sizeof(type);
			w_memcpy_s(&group, sizeof(group), &header_buf[i_start], sizeof(group));
			i_start += sizeof(group);
			w_memcpy_s(&run_id, sizeof(run_id), &header_buf[i_start], sizeof(run_id));
			i_start += sizeof(run_id);
			//w_memcpy_s(&desc, sizeof(desc), &header_buf[i_start], sizeof(desc));
			// This is done to remove possible system dependicies on whether char/uchar
			// is use to represent a standard char
			for (int i = 0; i < DESC_LEN; ++i)
			{
				if (!allowable_ascii_char(header_buf[i_start + 1]))
				{
					corrupt_desc = true;
					n = -2;
					return n;
				}
				else
				{
					desc[i] = header_buf[i_start + 1];
				}
			}
			i_start += sizeof(desc);
			desc[DESC_LEN - 1] = '\0';
			//get data
			data_len = buf_sz - i_start;
			data.resize(data_len, '\0');
			if (data_len > 0) {
				n = w_recvall(sockfd, &data[0], &data_len);
				if (data_len != buf_sz - header_sz)
				{
					n = -2;
					cerr << "NetPackage::recv error reading data: expected" << buf_sz - header_sz
						<< " bytes, but received " << data_len << "bytes" << endl;
				}
			}
		}
		if (n > 1) { n = 1; }
	}
	catch (exception& e)
	{
		cerr << "NetPackage::recv error" << endl;
		cerr << e.what() << endl;;
		n = -1;
	}
	catch (...)
	{
		n = -1;
	}
	return n;  // -2 on corrupt read, -1 on failure, 0 on a close connection or 1 on success
}

void NetPackage::print_header(std::ostream &fout)
{
	fout << "NetPackage: type = " << int(type) <<", group = " << group << ", run_id = " << run_id << ", description = " << desc <<
		", data package size = " << data.size() << endl;
}

//template std::string NetPackage::extract_string< std::vector<int8_t>::iterator>(std::vector<int8_t>::iterator first, std::vector<int8_t>::iterator last);
template std::vector<int8_t> NetPackage::pack_string< std::string::iterator>(std::string::iterator first, std::string::iterator last);
template std::vector<int8_t> NetPackage::pack_string< std::string::const_iterator>(std::string::const_iterator first, std::string::const_iterator last);
