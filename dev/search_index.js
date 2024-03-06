var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MyFirstPackage","category":"page"},{"location":"#MyFirstPackage","page":"Home","title":"MyFirstPackage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MyFirstPackage.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MyFirstPackage]","category":"page"},{"location":"#MyFirstPackage.AbstractLBConfig","page":"Home","title":"MyFirstPackage.AbstractLBConfig","text":"AbstractLBConfig{D, N}\n\nAn abstract type for lattice Boltzmann configurations.\n\nD: the dimension of the lattice\nN: the number of velocities\n\n\n\n\n\n","category":"type"},{"location":"#MyFirstPackage.D2Q9","page":"Home","title":"MyFirstPackage.D2Q9","text":"D2Q9 <: AbstractLBConfig{2, 9}\n\nA lattice Boltzmann configuration for 2D, 9-velocity model.\n\n\n\n\n\n","category":"type"},{"location":"#MyFirstPackage.curl-Union{Tuple{Array{Point2D{T}, 2}}, Tuple{T}} where T","page":"Home","title":"MyFirstPackage.curl","text":"curl(u::AbstractMatrix{Point2D{T}})\n\nCompute the curl of the momentum field in 2D, which is defined as:\n\nu_yxu_xy\n\n\n\n\n\n","category":"method"},{"location":"#MyFirstPackage.momentum-Tuple{MyFirstPackage.AbstractLBConfig, MyFirstPackage.Cell}","page":"Home","title":"MyFirstPackage.momentum","text":"momentum(lb::AbstractLBConfig, rho::Cell)\n\nCompute the momentum of the fluid from the density of the fluid.\n\n\n\n\n\n","category":"method"},{"location":"#MyFirstPackage.step!-Tuple{LatticeBoltzmann}","page":"Home","title":"MyFirstPackage.step!","text":"step!(lb::LatticeBoltzmann)\n\nPerform a single step of the lattice Boltzmann simulation.\n\n\n\n\n\n","category":"method"}]
}
