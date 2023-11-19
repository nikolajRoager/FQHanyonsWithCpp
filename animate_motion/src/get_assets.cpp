#include"get_assets.hpp"
#include<filesystem>

using namespace std;
namespace fs = std::filesystem;

//This better be the only place I use OS specific code, yes Microsoft provides a lot of nice c++ headers for Windows, and the same is true for Linux and Mac, but I refuse to sacrifice portability, not the least because that would limit me to Linux only.

#ifdef __unix__

fs::path get_assets()
{
	//I love Linux, but one thing it does not have is an application data folder, it really should have
	return fs::path("assets");
}
#elif _WIN32

fs::path get_assets()
{

	return fs::path("C:\\Program\ Files\ (x86)\\PGC_windows_export\\assets");
}

#endif

