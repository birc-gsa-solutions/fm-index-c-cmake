cmake_minimum_required (VERSION 3.21)
project (PROJECT)
enable_testing()

add_subdirectory(src)

add_executable(fm src/fm.c src/fasta.c src/fastq.c src/sam.c)
target_link_libraries(fm cstr)
