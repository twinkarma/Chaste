

# The auto-generated source files
set(config_dir ${Chaste_BINARY_DIR}/src)

#Generate Version.cpp
set(CHASTE_WC_MODIFIED "false")
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/ReleaseVersion.txt")
if(EXISTS "${version_file}")
    file(STRINGS "${version_file}" v_data)
    list(GET v_data 0 full_version)
    string(REPLACE "." ";" full_version_list "${full_version}")
    list(LENGTH full_version_list len)
    math(EXPR length ${len}-1)
    list(GET full_version_list ${length} chaste_revision)
    message("Chaste Release Full Version = ${full_version}, Revision = ${chaste_revision}")
else()
    # ReleaseVersion file not found, obtain revision information from SVN
	# The following requires a proper command-line svn client to be installed, not
	# just an ordinary shell extension like TortoiseSVN'
	# Install SlikSVN or the distribution from Collabnet (if you don't mind registering)
	find_package(Subversion)
	if(SUBVERSION_FOUND)
	   Subversion_WC_INFO("${Chaste_SOURCE_DIR}" chaste)
	   set(chaste_revision "${chaste_WC_REVISION}")
	   message("Current Chaste SVN Revision = ${chaste_WC_REVISION}. Chaste Last Changed Revision = ${chaste_WC_LAST_CHANGED_REV}")
	   if(${chaste_WC_REVISION} EQUAL ${chaste_WC_LAST_CHANGED_REV})
	      set(chaste_WC_MODIFIED "false")
	   else()
	      set(chaste_WC_MODIFIED "true")
	   endif()
	else(SUBVERSION_FOUND)
       # ReleaseVersion file not found, obtain revision information from Git
       find_package(Git REQUIRED)
       Git_WC_INFO("${CHASTE_SOURCE_ROOT}" chaste)
       set(chaste_revision "${chaste_WC_REVISION}")
       message("Current Chaste Git Revision = ${chaste_WC_REVISION}. Chaste Modified = ${chaste_WC_MODIFIED}")
    endif(SUBVERSION_FOUND)
endif()

#The generated timekeeper.cpp code below keeps track of build timestamp.
#It is built and executed prior to starting a build and prints the timestamp
#in a given format. This timestamp is used by Version.cpp, which is also auto-generated.
file(WRITE ${config_dir}/timekeeper.cpp
"#include <iostream>
#include <fstream>
#include <ctime>
int main( )
{
   time_t now = time(0);
   tm* loc_time = localtime(&now);
   char buffer[80];
   strftime(buffer, 80, \"%a, %d %b %Y %H:%M:%S +0000\", loc_time);

   std::ofstream timestampFile;
   timestampFile.open (\"build_timestamp\");
   timestampFile << buffer;
   timestampFile.close();
   return 0;
}
")

add_executable(timekeeper "${config_dir}/timekeeper.cpp")


#string(TIMESTAMP build_time)
set(build_time "set_in_configure_step")
execute_process(COMMAND ${XSD_EXECUTABLE} "--version" ERROR_VARIABLE xsd_version_full)
string(REGEX MATCH "^XML Schema Definition Compiler" xsd_version_2 "${xsd_version_full}")
string(REGEX MATCH "^CodeSynthesis XSD XML Schema to C\\+\\+ compiler" xsd_version_3 "${xsd_version_full}")
if (xsd_version_2) 
    set(xsd_version "2")
elseif (xsd_version_3)
    set(xsd_version "3")
else()
    set(xsd_version "undertermined")
endif()

set(time_size 80)
set(time_format "%a, %d %b %Y %H:%M:%S +0000")
set(project_versions "\"TODO\"")


find_package(PythonInterp)
execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" "from CheckForCopyrights import current_notice; print current_notice"
 WORKING_DIRECTORY "${Chaste_SOURCE_DIR}/python/infra"
 OUTPUT_VARIABLE licence)
string(REPLACE "\nThis file is part of Chaste.\n" "" licence "${licence}")
string(REPLACE "\n" "\\n" licence "${licence}")
set(quote "\"")
string(REPLACE ${quote} "\\${quote}" licence "${licence}")


set (Chaste_VERSION_MAJOR 3)
set (Chaste_VERSION_MINOR 3)



# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "global/src/Version_cmake.cpp.in"
  ${config_dir}/Version.cpp
  )


set(additional "")

if(MSVC)
    set(additional
    "
    //This has been put here to satisfy an MSVC linker issue 
    int __cdecl _purecall(void){return 0;}
    ")
endif()

configure_file (
  "global/src/ChasteBuildInfo_cmake.cpp.in"
  ${config_dir}/ChasteBuildInfo.cpp
  )

add_custom_target(generateTimestamp ALL
    COMMAND "$<TARGET_FILE:timekeeper>"
    COMMAND ${CMAKE_COMMAND} -P
        ${Chaste_SOURCE_DIR}/cmake/Modules/ChasteUpdateBuildTime.cmake
    WORKING_DIRECTORY "${config_dir}"
    DEPENDS timekeeper
    COMMENT "Generating Build Timestamp"
    VERBATIM
)


#add_custom_target(generateTimestamp ALL DEPENDS updateVersionCpp)

 

