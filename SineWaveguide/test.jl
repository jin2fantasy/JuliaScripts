function test(x::Type{Val{1}})
    println(x)
end

function test(x::Type{Val{2}})
    println(x)
end

c = 2
test(Val{c})

const aaa = rand(100,100);
@time aaa .== cc;

    bb = similar(aaa)
@time begin
    @inbounds for i in eachindex(aaa)
        bb[i] = (aa[i] == 0)
    end
end

@time map(.==, aa, zeros(100,100))
