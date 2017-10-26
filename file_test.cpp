#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>

int main(int argc, char ** argv)
{
    using namespace std;

    ofstream output_file("hello_file");

    if(!output_file.is_open())
    {
        cerr << "\tCannot open file hello_file" << endl;
    }
    else
    output_file << "Hello, world!" << endl;

    ifstream input_file("input_file");
    if(!input_file.is_open())
    cerr << "\tCannot open file: input_file" << endl;
    else
    while(!input_file.eof())
    {
        string str;
        input_file >> str;
        cout << str;
    }

    cout << "Command-line arguments:" << endl;
    for(int i = 0; i < argc; ++i)
    cout << argv[i] << endl;

    uint p = stoi(string(argv[1]), nullptr, 0);
    cout << "p = " << p << endl;

    vector<uint> vec(3);
    vec = {1, 2, 3};

    return 0;
}