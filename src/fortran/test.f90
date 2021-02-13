program test

print*, mod1(0,10)
print*, mod1(1,10)
print*, mod1(10,10)
print*, mod1(11,10)

contains

function mod1( x, y)
integer :: x
integer :: y

mod1 = modulo(x-1, y) + 1

end function mod1


end
