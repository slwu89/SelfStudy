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

//test making a struct
struct point {
    int x;
    int y;
    point(int xIn, int yIn) : x(xIn), y(yIn) {};
};

//main routine
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
    
    /*
     Pointers!!
     */
    //modifying what a pointer points to
    int i, j;
    i = 5; j = 2;
    int const *p = &i; //p points to i, which is a const int
    p = &j;  //p now points to j, which is not a const int
    
    //modifying a pointer itself
    int k = 5;
    int * const l = &k;
    std::cout << "memory address" << l << " " << *l << std::endl;
    *l = 42; //this is valid, now l points to the literal 42
    
    /*
     struct
     */
    point p1(4, 5);
    std::cout << "p1 info: " << p1.x << " " << p1.y << std::endl;
    
    return 0;
}
