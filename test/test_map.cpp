// #include <catch2/catch.hpp>
// #include <iostream>
// #include <typeinfo>
// #include <typeindex>
// using namespace Catch::Benchmark;


// TEST_CASE("map")
// {
// struct S {
// public:
//     int id;
//     std::string name;

//     S (int id, std::string name) : id(id), name(name) {}

//     bool operator<(const S& other) const {
//         if (id != other.id) {
//             return id < other.id;
//         }
//         return name < other.name;
//     }
// };

//     std::map<S, std::string> myMap;

//     S obj1(1, "AA");
//     S obj2(2, "BBB");
//     S obj3(1, "CCCCCC");


//     myMap[obj1] = "First";
//     myMap[obj2] = "Second";
//     myMap[obj3] = "Third";

//     for (const auto& entry : myMap) {
//         std::cout << "Key: (" << entry.first.id << ", " << entry.first.name << "), Value: " << entry.second << std::endl;
//     }



//     BENCHMARK("stdmap"){

//     };
// }

// TEST_CASE("unodered_map")
// {
// struct S  {
// public:
//     int id;
//     std::string name;

//     S(int id, std::string name) : id(id), name(name) {}

//        bool operator==(const S& other) const {
//         return id == other.id && name == other.name; // Выносится и заполняется через кортежи, tuple
//     }


// };
    


// struct Hash_fn {
//     std::size_t operator()(const S& obj) const {
//         std::size_t h1 = std::hash<int>()(obj.id);
//         std::size_t h2 = std::hash<std::string>()(obj.name);
//         return h1 ^ (h2>>5) ^ 2 ; // еще смещается и домнажается на константу.
//     }
// };
    
//     std::unordered_map<S, std::string, Hash_fn> myMap;

//     S obj1(1, "A");
//     S obj2(2, "B");
//     S obj3(1, "C");

//     std::type_index t(typeid(obj1));
//     // std::cout <<t.name();

//     myMap[obj1] = "First";
//     myMap[obj2] = "Second";
//     myMap[obj3] = "Third";

//     for (const auto& entry : myMap) {
//         std::cout << "Key: (" << entry.first.id << ", " << entry.first.name << "), Value: " << entry.second << std::endl;
//     }


//     BENCHMARK("stdunodered_mapmap"){

//     };
// }

