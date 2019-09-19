program test2

!#define TEST

#ifdef TEST
    print *, "This should print if TEST is defined"
#else 
    print *, "This should print if TEST is not defined"
#endif

end program test2


