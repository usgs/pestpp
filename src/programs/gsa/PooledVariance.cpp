#include "PooledVariance.h"
#include <regex>

using namespace std;

PooledVariance::PooledVariance()
{
}
void PooledVariance::process_pva_file(std::ifstream &fin)
{
	regex r("regex\"[^\"]+\"", regex_constants::icase);
	string line;
	while (getline(fin, line))
	{
		cout << line << endl;
		cmatch res_match;
		regex_search(line.c_str(), res_match, r);
		if (!res_match.empty())
		{
			for (auto x : res_match) std::cout << x << " ";
			cout << endl;
		}
	}
	cout << "################################" << endl;
//for (auto &i : res_match) cout << i << endl;
}

PooledVariance::~PooledVariance()
{
}
