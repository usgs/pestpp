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
#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

class FileManager
{
public:
	FileManager();
	FileManager(const std::string &_base_filename, const std::string &_directory="");
	void initialize_path(const std::string &_base_filename = "", const std::string &_directory = "");
	void open_default_files(bool restart = false);
	std::string build_filename(const std::string &ext, bool include_dir = false);
	//std::string get_full_filename(const std::string &tag);
	std::string get_base_filename() {return pest_base_filename;}
	std::ofstream &rec_ofstream();
	std::ofstream &sen_ofstream();
	std::ofstream &open_ofile_ext(const std::string &extension, std::ios_base::openmode mode = std::ofstream::out);
	std::ofstream &open_ofile_local(const std::string &tag, const std::string &filename, std::ofstream::openmode mode = std::ofstream::out);
	std::ofstream &open_ofile_absolute(const std::string &tag, const std::string &filename, std::ofstream::openmode mode = std::ofstream::out);
	std::ifstream &open_ifile_ext(const std::string &extension, std::ifstream::openmode mode = std::ifstream::in);
	std::ifstream &open_ifile_local(const std::string &tag, const std::string &filename, std::ifstream::openmode mode = std::ifstream::in);
	std::ifstream &open_ifile_absolute(const std::string &tag, const std::string &filename, std::ifstream::openmode mode = std::ifstream::in);
	std::fstream &open_iofile_ext(const std::string &extension, std::ios_base::openmode mode = std::fstream::in | std::fstream::out);
	std::fstream &open_iofile_local(const std::string &tag, const std::string &filename, std::fstream::openmode mode = std::fstream::in | std::fstream::out);
	std::fstream &open_iofile_absolute(const std::string &tag, const std::string &filename, std::fstream::openmode mode = std::fstream::in | std::fstream::out);
	void close_file(const std::string &extension);
	std::ofstream &get_ofstream(const std::string &tag);
	std::ifstream &get_ifstream(const std::string &tag);
	std::fstream &get_fstream(const std::string &tag);
	~FileManager(void);
private:
	std::string directory;
	std::string pest_base_filename;
	std::map<std::string, std::ofstream*> ofile_map;
	std::map<std::string, std::ifstream*> ifile_map;
	std::map<std::string, std::fstream*> iofile_map;
	std::map<std::string, std::string> filename_map;
};

#endif /* FILEMANAGER_H_ */
