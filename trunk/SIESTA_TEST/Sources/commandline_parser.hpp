//------------------------------------------------------------------------------
//  The @header2, @begin_table, @item3 and @end_table commands are custom
//  defined commands in Doxygen.in. They are defined under ALIASES. For the page
//  created here, the 80 column limit is exceeded. Arguments of aliases are
//  separated by ','. If you intended ',' to be a string you must use an escaped
//  comma '\,'.
//
///  @page siesta_test_cl_parsing_sec Command Line Arguments
///
///  @tableofcontents
///
///  @section siesta_test_cl_parsing_intro Introduction
///  This contains a description of the command line arguments. All arguments
///  take the form of
///
///  @fixed_width{-arg=value}
///
///  @section siesta_test_cl_parsing_arg_sec Command Line Arguments
///  @header2{Argument, Takes Value, Discription}
///  @begin_table
///     @item3{@fixed_width{-h},            N, Displays the help text and exits the program.}
///     @item3{@fixed_width{-wout_file},    Y, Specify the mgrid input file name.}
///     @item3{@fixed_width{-restart_file}, Y, Specify the wout input file name.}
///  @end_table
///
///
///  @section siesta_test_cl_pasring_prog_ref_sec Programmers Reference
///  Reference material for the coding to implement command line parsing is found
///  in the @ref commandline_parser class.
//------------------------------------------------------------------------------
//******************************************************************************
///  @file commandline_parser.hpp
///  @brief Contains classes to interpolate full and half grid quanities.
//******************************************************************************

#ifndef commandline_parser_hpp
#define commandline_parser_hpp

#include <iostream>
#include <map>
#include <string>

///  Type for the argument map.
typedef std::map<std::string, std::string> arg_map;
///  Type for a key value pair.
typedef std::pair<std::string, std::string> arg_element;

//------------------------------------------------------------------------------
///  @brief A radial quantity.
//------------------------------------------------------------------------------
class commandline_parser {
public:
/// Parsed commands.
    const arg_map commands;

//------------------------------------------------------------------------------
///  @brief Factory method to parse the commandline and produce the arguments.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Commandline strings.
///  @returns A constructed map of commandline argument key value pairs.
//------------------------------------------------------------------------------
    static arg_map parse_commands(const size_t argc, const char * argv[]) {
        if (argc == 0) {
            help();
        }

        arg_map commands;


        for (size_t i = 1; i < argc; i++) {
            std::string arg(argv[i]);

            if (arg == "-h") {
                help();
            } else {
                size_t eqpos = arg.find('=');
                if (eqpos != std::string::npos) {
                    std::string key = arg.substr(0, eqpos);
                    std::string value = arg.substr(eqpos + 1, std::string::npos);

                    commands.insert(arg_element(key, value));
                } else {
                    commands.insert(arg_element(arg, ""));
                }
            }
        }

        return commands;
    }

//------------------------------------------------------------------------------
///  @brief Construct a commandline_parser object by pasring command line
///  arguments.
///
///  @param[in] argc Number of commandline arguments.
///  @param[in] argv Commandline strings.
//------------------------------------------------------------------------------
    commandline_parser(const size_t argc, const char * argv[]) :
    commands(parse_commands(argc, argv)) {}

//------------------------------------------------------------------------------
///  @brief Check if command arg was set.
///
///  @param[in] key Commandline key to check.
///  @returns True if the key was set.
//------------------------------------------------------------------------------
    bool is_set(const std::string &key) const {
        return commands.find(key) != commands.end();
    }

//------------------------------------------------------------------------------
///  @brief Check if command arg was set.
///
///  @param[in] key Commandline key to check.
///  @returns True if the key was set.
//------------------------------------------------------------------------------
    std::string get(const std::string &key) const {
        if (!this->is_set(key)) {
            help();
        }

        return commands.at(key);
    }

//------------------------------------------------------------------------------
///  @brief Display help.
//------------------------------------------------------------------------------
    static void help() {
//                   "                                        ''                                      "
        std::cout << "                                                                                " << std::endl;
        std::cout << "                                    SIESTA TEST                                 " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "Usage: xsiesta_test [-arg][=option] ...                                         " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "Options:                                                                        " << std::endl;
        std::cout << "All options are displayes as [arg][takesoption][Discription]                    " << std::endl;
        std::cout << "  -h            N Display this information.                                     " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -wout_file    Y Path to the VMEC wout file.                                   " << std::endl;
        std::cout << "                                                                                " << std::endl;
        std::cout << "  -restart_file Y Path to the SIESTA restart file.                              " << std::endl;
        std::cout << std::endl;

        exit(0);
    }
};

#endif
