/**
   @file grid_test.cpp

   Tests generating function grids
*/

#include "sput.h"
#include <sstream>

#include <grid.hpp>

struct SomeFunction {
    virtual ~SomeFunction() {}
    virtual double operator()(double x) const =0;
};

struct Square : public SomeFunction {
    virtual double operator()(double x) const { return x*x; }
};

struct Cube : public SomeFunction {
    virtual double operator()(double x) const { return x*x*x; }
};
    

typedef Grid<SomeFunction> MyGrid;

static void sanity()
{
    {
        Square s;
        const SomeFunction& f=s;
        sput_fail_unless(f(2.)==4., "Square function sanity check");
    }
    {
        Cube c;
        const SomeFunction& f=c;
        sput_fail_unless(f(2.)==8., "Cube function sanity check");
    }
}

static void ctor()
{
    MyGrid grid;
    const MyGrid& cgrid=grid;
    
    sput_fail_unless(cgrid.columns()==0, "Initial columns");
    sput_fail_unless(cgrid.xmin()==0, "Initial xmin");
    sput_fail_unless(cgrid.xmax()==0, "Initial xmax");
    sput_fail_unless(cgrid.rows()==512, "Initial rows");

    Square sf;
    Cube cf;
    grid.set_xmin(-5)
        .set_xmax(5)
        .set_rows(1024)
        .add_column(&sf)
        .add_column(&cf);

    sput_fail_unless(cgrid.columns()==2, "Set columns");
    sput_fail_unless(cgrid.xmin()==-5.0, "Set xmin");
    sput_fail_unless(cgrid.xmax()==5.0, "Set xmax");
    sput_fail_unless(cgrid.rows()==1024, "Set rows");
}

static void values()
{
    Square sf;
    Cube cf;
    MyGrid grid;

    const double expected[][3]={
        { -2, 4, -8 },
        {  0, 0,  0 },
        {  2, 4,  8 }
    };
    const unsigned int nrows=sizeof(expected)/sizeof(*expected);
    const unsigned int ncols=sizeof(expected[0])/sizeof(*(expected[0]));
    
    grid.set_xmin(-2)
        .set_xmax(2)
        .set_rows(nrows)
        .add_column(&sf)
        .add_column(&cf);

    std::ostringstream out;
    out << grid;

    std::istringstream ins(out.str());
    std::string line;
    for (unsigned int i=0; i<nrows; ++i) {
        getline(ins, line);
        sput_fail_if(!ins, "can read line");
        if (!ins) return;
        std::istringstream line_ins(line);
        for (unsigned int j=0; j<ncols; ++j) {
            double x;
            line_ins >> x;
            sput_fail_if(!line_ins, "can read");
            if (!line_ins) return;
            sput_fail_unless(expected[i][j]==x, "expected value");
        }
        sput_fail_unless(line_ins.eof(), "No extra characters in line");
    }
    getline(ins, line);
    sput_fail_unless(!ins, "No extra lines in output");
}

int main(int argc, char **argv)
{
    sput_start_testing();
	
    sput_enter_suite("test_grid_test");

    sput_run_test(sanity);
    sput_run_test(ctor);
    sput_run_test(values);

    sput_finish_testing();
    return sput_get_return_value();
}
