#include <exception>
#include <string>

namespace error {

	class input_err : public std::exception {
		std::string err_msg;

		public:
		input_err(std::string msg) {
			err_msg = "Excpected A, C, G, T but received '" + msg + "'.";

		};
		
		std::string describe();
	};

}
