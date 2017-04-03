//Thoth Gunter
//

void parse_argv( std::string& tmp_options, std::vector<std::string>& options )
{
  bool more_to_parse = true;
  int index = 0;
  int iter = 0;
  std::cout << tmp_option << std::endl;
	while(more_to_parse == true){
		std::string token = tmp_option.substr(index, tmp_option.find(" ", index, 1) - index );
		options.push_back(token);
		iter += 1;
		if( tmp_option.find(" ", index, 1) == std::string::npos ) more_to_parse = false;
		i = tmp_option.find(" ", index, 1) + 1 ;
	}
	
};
