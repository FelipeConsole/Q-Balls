
using Pkg
Pkg.add.(["CSV","DataFrames","PlutoUI","PyPlot"])

using CSV

using DataFrames

using PlutoUI

Pkg.add("PyPlot")
using PyPlot


begin
    csv_data = CSV.File("phi0_1.csv")
    data = DataFrame(csv_data)
end


begin
    csv_data2 = CSV.File("phi0_1_v2.csv")
    data2 = DataFrame(csv_data2)
end


begin
    csv_data3 = CSV.File("phi_0=0_4.csv")
    data4 = DataFrame(csv_data3)
end

Pkg.add("DataFrames")
