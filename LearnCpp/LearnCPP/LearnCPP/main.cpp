//
//  main.cpp
//  LearnCPP
//
//  Created by Sean Wu on 11/17/16.
//  Copyright © 2016 Sean Wu. All rights reserved.
//

#include <iostream>

//test making a typedef
typedef char testTypeDef[3];

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    
    //what the hell are enumerators?
    enum color {red,blue,green,black};
    color x;
    x = red;
    color x1 = blue;
    color x2 = green;
    std::cout << x << " " << x1 << " " << x2 << std::endl;
    
    //test out the typedef and print first element
    testTypeDef testTypeDefVar[] = {"a","b","c"};
    std::cout << testTypeDefVar[0] << std::endl;
    
    
    return 0;
}
