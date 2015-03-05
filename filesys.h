#include <iostream>
#include <list>
#include <string>

#ifndef NOSTLFS

#include <filesystem>

using namespace std;
using namespace std::tr2::sys;

namespace xmv
{
	//检索指定目录中的所有文件
	string GetFileName(string path);

	template<class ListOfString>
	void GetFiles(string folderName, ListOfString & filepathList, string extName = "*.*",  bool findSubFolder = false)
	{
		path folerPath(folderName);

		if ( ! exists(folerPath) ) return; 

		bool noFilter = extName == "*" || extName == "*.*" || extName == "" || extName.empty();
		if (! noFilter)
		{
			if (extName.size() >= 2 && extName.substr(0, 2) == "*.") extName = extName.substr(1);
			if (extName.size() >= 1 && extName.substr(0, 1) != ".") extName = "." + extName;
		}

		directory_iterator end_itr; // 缺省构造生成一个结束迭代器 
		for (directory_iterator itr(folerPath); itr != end_itr; ++itr) 
		{
			string fn = itr->path();
			if (is_directory(itr->status())) 
			{ 
				if (findSubFolder)
				{
					GetFiles(fn, filepathList, extName, findSubFolder); //递归查找 
				}
			} 
			else
			{ 
				if (noFilter)
					filepathList.push_back(fn);
				else
				{
					if (itr->path().extension() == extName)
						filepathList.push_back(fn);
				}
			}
		}
	}

#undef CreateDirectory
	void CreateDirectory(string directoryName);
	string StartupPath();
}


#endif