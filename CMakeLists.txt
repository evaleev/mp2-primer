project(mp2-primer)

add_executable(mp2 main.cc mp2.cc hartree-fock.cc)
target_link_libraries(mp2 int2)
