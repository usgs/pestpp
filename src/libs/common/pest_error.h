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

#ifndef PEST_ERROR_H_
#define PEST_ERROR_H_
#include <stdexcept>
#include <string>

using namespace std;

class PestError : public std::runtime_error {
public:
	PestError(const string &_message="") : runtime_error(_message), message(_message) {};
	virtual ~PestError() throw (){}
	void add_front(const string &s) {message = s + message;}
	void add_back(const string &s) {message += s;}
	void raise() {throw *this;}
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
protected:
	string message;
};

class PestConversionError : public PestError {
public:
	PestConversionError(const string &_str, const string &_message="")
		: PestError(_message) , str(_str){
		message = string("PestConversionError:  Error processing: \"") + _str + "\"" + message;
	}
	virtual ~PestConversionError() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string str;
};

class PestFileError : public PestError {
public:
	PestFileError(const string &_filename, const string &_message="")
		: PestError(_message) , filename(_filename){
		message = string("PestFileError:  Error opening file: \"") + filename + "\"" + message;
	}
	virtual ~PestFileError() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string filename;
};

class PestFileErrorAccess : public PestError {
public:
	PestFileErrorAccess(const string &_filename, const string &_message="")
		: PestError(_message) , filename(_filename){
		message = string("PestFileError:  Error accessing file: \"") + filename + "\"" + message;
	}
	virtual ~PestFileErrorAccess() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string filename;
};

class PestParsingError : public PestError {
public:
	PestParsingError(const string &_line, const string &_message="")
		: PestError(_message) , line(_line){
		message = string("PestParsingError:  Error parsing: \"") + line + "\"  " + message;
	}
	virtual ~PestParsingError() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string line;
};


class PestIndexError : public PestError {
public:
	PestIndexError(const string &_str, const string &_message="")
		: PestError(_message) , str(_str){
		message = string("PestIndexError:  Invalid index: \"") + _str + "\"" + message;
	}
	virtual ~PestIndexError() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string str;
};

class PestCommandlineError : public PestError {
public:
	PestCommandlineError(const string &_str, const string &_message="")
		: PestError(_message) , str(_str){
		message = string("Invalid command line: \"") + _str + "\"" + message;
	}
	virtual ~PestCommandlineError() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string str;
};


#endif /* PEST_ERROR_H_ */
