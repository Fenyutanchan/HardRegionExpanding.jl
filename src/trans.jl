using JLD2
using SymEngine
# A1=load("amp1.jld2","amp_lorentz_list")
# A0=load("amp1.jld2","amp_color_list")
# A2=load("amp1.jld2","kin_relation")
# A3=load("amp1.jld2","loop_den_list")
# A4=load("amp1.jld2","signed_symmetry_factor")
# amp=load("amp1.jld2")



# A10=map(Basic,amp["amp_color_list"])
# A11=map(Basic,amp["amp_lorentz_list"])
# A12=map(Basic,amp["loop_den_list"])
# A13=Basic(amp["signed_symmetry_factor"])



# amp_expr=sum(A10.*A11*prod(A12))*A13


# for rel in A2 
#     println(rel, rel[1])
# end

# amp_expr1=amp_expr

# for rel in A2 
#     amp_expr1=subs(amp_expr1,Basic(rel[2])=>Basic(rel[1]))
# end

process= Basic( replace(String(chop(filter(endswith(".yaml"),readdir(".."))[1],tail=5)),"_"=>""))
diagramnum=length(filter(endswith(".jld2"),readdir()))
println(process,"gogogo")

# println(amp_expr)
function ampcombine(amp_file::String,diagramnum::Int64)
    amp=load(amp_file)
    A10=map(Basic,amp["amp_color_list"])
    A11=map(Basic,amp["amp_lorentz_list"])
    A12=map(Basic,amp["loop_den_list"])
    A13=Basic(amp["signed_symmetry_factor"])
    A14=amp["kin_relation"]
    amp_expr=sum(A10.*A11*prod(A12))*A13
    for rel in A14
        amp_expr=subs(amp_expr,Basic(rel[2])=>Basic(rel[1]))
    end
    return amp_expr
end

#ampcombine("amp1.jld2")
form_list_file= open( "$process.frm", "w" )
write(form_list_file," #-\n#include contractor.frm\n#include hardregion.frm\n#define diagramnum \"$diagramnum\"\n")
for index in 1:diagramnum
   sus=ampcombine("amp$index.jld2",diagramnum)
   write( form_list_file, "Local expr$index=$sus;\n")
   write(form_list_file,"#call hardregion(1,1)\n #call ArrangeTrace()\n")
end
write(form_list_file,"L amp$process=expr1+...+expr`diagramnum';\n print;\n .end")
close(form_list_file)