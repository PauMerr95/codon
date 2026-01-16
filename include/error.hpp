#include <exception>
#include <string>

namespace err {

	class input_err : public std::exception {
		std::string err_msg;

		input_err(std::string msg);
		
		std::string describe();
	};
}
