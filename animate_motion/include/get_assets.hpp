#pragma once

#include<filesystem>

using namespace std;
namespace fs = std::filesystem;


//Calls a OS specific function to get the assets directory.
fs::path get_assets();//Note that I am intentionally shielding any OS specific code in just the one file get_assets.cpp for portability
