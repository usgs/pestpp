
/*


	This file is part of PEST++.

	PEST++ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	PEST++ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/
#include "FileManager.h"
#include "utilities.h"
#include "OutputFileWriter.h"
#include "system_variables.h"

using namespace pest_utils;

FileManager::FileManager()
	: directory(""), pest_base_filename("")
{
}
FileManager::FileManager(const string &_base_filename, const string &_directory)
  : directory(strip_cp(_directory)), pest_base_filename(strip_cp(_base_filename))
{
	initialize_path(_base_filename, _directory);
}

void FileManager::initialize_path(const string &_base_filename, const string &_directory)
{
	directory = strip_cp(_directory);
	pest_base_filename = strip_cp(_base_filename);
}

void  FileManager::open_default_files(bool restart)
{
	close_file("rec");
	if (restart)
	{
		open_ofile_ext("rec", fstream::app);
		open_ofile_ext("rst", fstream::app);
	}
	else
	{
		open_ofile_ext("rec");
		open_ofile_ext("rst");
	}


}



FileManager::~FileManager(void)
{
	for (auto &i : ofile_map)
	{
		if(i.second->is_open()) i.second->close();
		delete i.second;
	}
	for (auto &i : ifile_map)
	{
		if(i.second->is_open()) i.second->close();
		delete i.second;
	}
	for (auto &i : iofile_map)
	{
		if(i.second->is_open()) i.second->close();
		delete i.second;
	}
}

string FileManager::build_filename(const string &ext, bool include_dir)
{
	if (!include_dir)
	{
		return  pest_base_filename + "." + strip_cp(ext);
	}
	else
	{
		return directory + OperSys::DIR_SEP + pest_base_filename + "." + strip_cp(ext);
	}
}

//string FileManager::get_full_filename(const string &tag)
//{
//	string name = "NOT_FOUND";
//	auto i = filename_map.find(tag);
//	if (i != filename_map.end())
//	{
//		name = i->second;
//	}
//	return directory + OperSys::DIR_SEP + name;
//}

ofstream &FileManager::open_ofile_absolute(const string &tag, const string &filename, ofstream::openmode mode)
{
	pair<map<string ,ofstream*>::iterator, bool> ret;
	ret = ofile_map.insert(pair<string, ofstream*>(tag, new ofstream()));
	ofstream &f_new = *(ret.first->second);
	if (ret.second != false  && !f_new.is_open())
	{
		f_new.open(filename, mode);
		filename_map.insert(pair<string, string>(tag, filename));
	}
	if (f_new.bad())
	{
		throw PestFileError(filename);
	}
	return f_new;
}

ofstream &FileManager::open_ofile_ext(const string &extension, ofstream::openmode mode)
{
		string filename = build_filename(extension);
		ofstream &f_new = open_ofile_absolute(extension, filename, mode);
		return f_new;
}

ofstream &FileManager::open_ofile_local(const string &tag, const string &filename, ofstream::openmode mode)
{
	string comp_filename = directory + OperSys::DIR_SEP + filename;
	ofstream &f_new = open_ofile_absolute(tag, filename, mode);
	return f_new;
}

ifstream &FileManager::open_ifile_absolute(const string &tag, const string &filename, ifstream::openmode mode)
{
	pair<map<string ,ifstream*>::iterator, bool> ret;
	ret = ifile_map.insert(pair<string, ifstream*>(tag, new ifstream()));
	ifstream &f_new = *(ret.first->second);
	if (ret.second != false  && !f_new.is_open())
	{
		f_new.open(filename, mode);
		filename_map.insert(pair<string, string>(tag, filename));
	}
	if (f_new.bad())
	{
		PestFileError pest_err(filename);
		cerr << pest_err.what();
		throw pest_err;
	}
	return f_new;
}

ifstream &FileManager::open_ifile_ext(const string &extension, ifstream::openmode mode)
{
		string filename = build_filename(extension);
		ifstream &f_new = open_ifile_absolute(extension, filename, mode);
		return f_new;
}

ifstream &FileManager::open_ifile_local(const string &tag, const string &filename, ifstream::openmode mode)
{
	string comp_filename = directory + OperSys::DIR_SEP + filename;
	ifstream &f_new = open_ifile_absolute(tag, filename, mode);
	return f_new;
}

fstream &FileManager::open_iofile_absolute(const string &tag, const string &filename, fstream::openmode mode)
{
	pair<map<string ,fstream*>::iterator, bool> ret;
	ret = iofile_map.insert(pair<string, fstream*>(tag, new fstream()));
	fstream &f_new = *(ret.first->second);
	if (ret.second != false  && !f_new.is_open())
	{
		// a file needs to exist before it can be opened it with read and write
		// permission.   So open it with write permission to crteate it, close
		// and then reopen it with read and write permission.
		f_new.open(filename, ios_base::out | ios_base::binary);
		f_new.close();
		f_new.open(filename, mode);
		filename_map.insert(pair<string, string>(tag, filename));
	}
	if (!f_new.good())
	{
		throw PestFileError(filename);
	}
	return f_new;
}

fstream &FileManager::open_iofile_ext(const string &extension, fstream::openmode mode)
{
		string filename = build_filename(extension);
		fstream &f_new = open_iofile_absolute(extension, filename, mode);
		return f_new;
}

fstream &FileManager::open_iofile_local(const string &tag, const string &filename, fstream::openmode mode)
{
	string comp_filename = directory + OperSys::DIR_SEP + filename;
	fstream &f_new = open_iofile_absolute(tag, filename, mode);
	return f_new;
}


void FileManager::close_all_files_containing(const std::string& partial_name)
{
    string lower = lower_cp(partial_name);
    vector<string> to_close;
    for (auto& fm : filename_map)
    {
        string lower_fm = lower_cp(fm.second);
        if (lower_fm.find(lower) != string::npos)
        {
            to_close.push_back(fm.first);
        }
    }
    for (auto tc : to_close)
        close_file(tc);
}

void FileManager::close_all_files_ending_with(const std::string& partial_name)
{
	string lower = lower_cp(partial_name);
	vector<string> to_close;
	for (auto& fm : filename_map)
	{
		//if (lower_cp(fm.second).find(lower) != string::npos)
		string lower_fm = lower_cp(fm.second);
		string t = lower_fm.substr(lower_fm.size() - partial_name.size(), lower_fm.size());
		if (t == lower)
		{
			to_close.push_back(fm.first);
		}
	}
	for (auto tc : to_close)
		close_file(tc);
}

void FileManager::close_file(const string &tag)
{
	{
		map<string, ofstream*>::iterator it;
		it = ofile_map.find(tag);
		if (it != ofile_map.end() && it->second->is_open())
		{
			it->second->close();
			delete it->second;
			ofile_map.erase(it);
		}
	}
	{
		map<string, ifstream*>::iterator it;
		it = ifile_map.find(tag);
		if (it != ifile_map.end() && it->second->is_open())
		{
			it->second->close();
			delete it->second;
			ifile_map.erase(it);
		}
	}
	{
		map<string, fstream*>::iterator it;
		it = iofile_map.find(tag);
		if (it != iofile_map.end() && it->second->is_open())
		{
			it->second->close();
			delete it->second;
			iofile_map.erase(it);

		}
	}
	filename_map.erase(tag);
}

ofstream &FileManager::get_ofstream(const string &tag)
{
	map<string, ofstream*>::iterator it;
	it = ofile_map.find(tag);
	if (it == ofile_map.end())
	{
		string filename = build_filename(tag);
		throw PestFileErrorAccess(filename);
	}
	return *(it->second);
}

fstream &FileManager::get_fstream(const string &tag)
{
	map<string, fstream*>::iterator it;
	it = iofile_map.find(tag);
	if (it == iofile_map.end())
	{
		string filename = build_filename(tag);
		throw PestFileErrorAccess(filename);
	}
	return *(it->second);
}

ifstream &FileManager::get_ifstream(const string &tag)
{
	map<string, ifstream*>::iterator it;
	it = ifile_map.find(tag);
	if (it == ifile_map.end())
	{
		string filename = build_filename(tag);
		throw PestFileErrorAccess(filename);
	}
	return *(it->second);
}


ofstream &FileManager::rec_ofstream()
{
	return get_ofstream("rec");
}

ofstream &FileManager::sen_ofstream()
{
	return get_ofstream("sen");
}
