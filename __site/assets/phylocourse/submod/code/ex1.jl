# This file was generated, do not modify it. # hide
function simulate(site, n, p)
    print(site)
    for i=1:n
        if rand() < p
            site = rand(setdiff("ATCG", site))
        end
        print(" --> ", site)
    end
end