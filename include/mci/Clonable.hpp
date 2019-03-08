#ifndef MCI_CLONABLE_HPP
#define MCI_CLONABLE_HPP

#include <memory>

// Simple template base for safely clonable interfaces
//
// Derive from this like
//
// class MyClass: Clonable<MyClass> {...}
//
// and implement method _clone which should return a raw pointer of type MyClass,
// pointing to a cloned version of the object. In effect, your class will have
// a safe public clone() method returning a unique pointer of type MyClass.
//
// Note: If you are implementing a child class of a clonable base, like
//
// class MyChild: public MyBase {...}
// MyBase: Clonable<MyBase>
//
// you might also want add another non-overriding public clone() method returning
// std::unique<MyChild>, instead of std::unique<MyBase>. This means calling clone()
// on a pointer of type MyBase returns unique MyBase pointer and clone() on MyChild
// returns unique MyChild pointer.

namespace mci
{
    template<typename T>
    struct Clonable
    {
    protected:
        virtual T * _clone() const = 0;

    public:
        virtual ~Clonable() = default;

        std::unique_ptr<T> clone() const
        {
            return std::unique_ptr<T>(_clone());
        }
    };
}  // namespace mci
#endif
