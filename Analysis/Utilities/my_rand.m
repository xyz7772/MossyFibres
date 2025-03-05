function [z] = my_rand(size, a, b)
    z = a + (b-a) * rand(size);
end