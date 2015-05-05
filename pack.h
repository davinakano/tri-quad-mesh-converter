#ifndef PACK_H
#define PACK_H

template<class T1, class T2>
class pack2
{
public:

    pack2()
    {

    }

    pack2(const pack2<T1, T2> &ref)
    {
        this->first = ref.first;
        this->second = ref.second;
    }

    pack2(const T1 &_first, const T2 &_second)
    {
        first = _first;
        second = _second;
    }

    ~pack2()
    {

    }

    bool operator<(const pack2<T1, T2> &right) const
    {
        return first < right.first;
    }

    bool operator>(const pack2<T1, T2> &right) const
    {
        return first > right.first;
    }

    bool operator==(const pack2<T1, T2> &right) const
    {
        return first == right.first;
    }

    bool operator!=(const pack2<T1, T2> &right) const
    {
        return first != right.first;
    }

    T1 first;
    T2 second;
};

template<class T1, class T2, class T3>
class pack3
{
public:

    pack3()
    {

    }

    pack3(const pack3<T1, T2, T3> &ref)
    {
        this->first = ref.first;
        this->second = ref.second;
        this->third = ref.third;
    }

    pack3(const T1 &_first, const T2 &_second, const T3 &_third)
    {
        first = _first;
        second = _second;
        third = _third;
    }

    ~pack3()
    {

    }

    bool operator<(const pack3<T1, T2, T3> &right) const
    {
        return first < right.first;
    }

    bool operator>(const pack3<T1, T2, T3> &right) const
    {
        return first > right.first;
    }

    bool operator==(const pack3<T1, T2, T3> &right) const
    {
        return first == right.first;
    }

    T1 first;
    T2 second;
    T3 third;
};

template<class T1, class T2, class T3, class T4>
class pack4
{
public:

    pack4()
    {

    }

    pack4(const pack4<T1, T2, T3, T4> &ref)
    {
        this->first = ref.first;
        this->second = ref.second;
        this->third = ref.third;
        this->fourth = ref.fourth;
    }

    pack4(const T1 &_first, const T2 &_second, const T3 &_third, const T4 &_fourth)
    {
        first = _first;
        second = _second;
        third = _third;
        fourth = _fourth;
    }

    ~pack4()
    {

    }

    bool operator<(const pack4<T1, T2, T3, T4> &right) const
    {
        return first < right.first;
    }

    bool operator>(const pack4<T1, T2, T3, T4> &right) const
    {
        return first > right.first;
    }

    bool operator==(const pack4<T1, T2, T3, T4> &right) const
    {
        return first == right.first;
    }

    bool operator!=(const pack4<T1, T2, T3, T4> &right) const
    {
        return first != right.first;
    }

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
};


typedef pack2<int, int> pack2int;
typedef pack3<int, int, int> pack3int;
typedef pack2<double, double> pack2double;

#endif // PACK_H
