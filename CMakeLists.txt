project(mp2-primer)

add_executable(hf hartree-fock.cc)
target_link_libraries(hf int2)
