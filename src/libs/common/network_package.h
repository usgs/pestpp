#ifndef NET_PACKAGE_H_
#define NET_PACKAGE_H_
#include <string>
#include <cstdint>
#include <vector>
#include <memory>

class NetPackage
{
public:
	static bool allowable_ascii_char(int8_t value);
	static bool check_string(const int8_t *data_src, size_t _size);
	static bool check_string(const std::vector<int8_t> &data_src, size_t index1, size_t _size);
	static std::string extract_string(const int8_t *data_src, size_t _size);
	static std::string extract_string(const std::vector<int8_t> &data_src, size_t index1, size_t _size);
	//template<class InputIterator>
	//static std::string extract_string(InputIterator first, InputIterator last);
	template<class InputIterator>
	static std::vector<int8_t> pack_string(InputIterator first, InputIterator last);
	//dont forget to add to pack_strings also!
	enum class PackType :uint32_t {
		UNKN, OK, CONFIRM_OK, READY, REQ_RUNDIR, RUNDIR, REQ_LINPACK, LINPACK, PAR_NAMES, OBS_NAMES,
		START_RUN, RUN_FINISHED, RUN_FAILED, RUN_KILLED, TERMINATE,PING,REQ_KILL,IO_ERROR,CORRUPT_MESG,
		DEBUG_LOOP,DEBUG_FAIL_FREEZE,START_FILE_WRKR2MSTR,CONT_FILE_WRKR2MSTR, FINISH_FILE_WRKR2MSTR,
		// APPEND ONLY: the value is serialized as uint32_t, so inserting anywhere above
		// renumbers every later type and two versions stop understanding each other.
		REQ_PARTIAL, PARTIAL_OBS};
	
	static int get_new_group_id();
	NetPackage(PackType _type=PackType::UNKN, int _group=-1, int _run_id=-1, const std::string &desc_str="");
	~NetPackage(){}

	std::vector<std::string> pack_strings;
	const static int DESC_LEN = 1001;
	/// Announced by an agent in the free-form text of its LINPACK and READY messages, and
	/// looked for by the master, to decide whether REQ_PARTIAL is safe to send to it. See
	/// AgentInfoRec::get_supports_partial().
	static constexpr const char* PARTIAL_CAPABILITY_TAG = "partial=1";
	/// Does this message text advertise partial-results support?
	///
	/// Both places the master asks live in the same if/else chain, and the answer has to be
	/// the same in both: an agent that says so during the handshake and an agent that says so
	/// after a run are equally able to answer REQ_PARTIAL. Sending REQ_PARTIAL to one that
	/// cannot is not a no-op - it treats the message as corrupt and kills the run - so this
	/// says no for anything it does not positively recognise, empty text included.
	static bool advertises_partial(const std::string& info_txt)
	{
		return info_txt.find(PARTIAL_CAPABILITY_TAG) != std::string::npos;
	}
	const static int NULL_DA_CYCLE = -9999;
	const static int FILE_TRANS_BUF_SIZE = 102400;
	std::pair<int,std::string> send(int sockfd, const void *data, int64_t data_len_l);
	std::pair<int,std::string> recv(int sockfd);
	void reset(PackType _type, int _group, int _run_id, const std::string &_desc);
	PackType get_type() const {return type;}
	int64_t get_run_id() const { return run_id; }
	int64_t get_group_id() const { return group; }
	std::string get_info_txt();
	const std::vector<int8_t> &get_data(){ return data; }
	void print_header(std::ostream &fout);


private:
	bool verbose;
	int64_t data_len;
	static int64_t last_group_id;
	PackType type;
	int64_t group;
	int64_t run_id;
	int8_t desc[DESC_LEN];
	static int8_t security_code[5];
	std::vector<int8_t> data;
};


#endif /* NET_PACKAGE_H_ */
