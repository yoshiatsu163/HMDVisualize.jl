const _Bond_Quality = 8
const _Single_Radius = 1.0
const _Double_Radius = 1/3
const _Triple_Radius = 1/5

const Single_Bond = Makie._mantle(Point3f(zeros(3)), Point3f(0,0,1), _Single_Radius, _Single_Radius, _Bond_Quality)

const Double_Bond = begin
    m1 = Makie._mantle(Point3f(-2/3, 0.0, 0.0), Point3f(-2/3, 0.0, 1.0), _Double_Radius, _Double_Radius, _Bond_Quality)
    m2 = Makie._mantle(Point3f( 2/3, 0.0, 0.0), Point3f( 2/3, 0.0, 1.0), _Double_Radius, _Double_Radius, _Bond_Quality)
    merge([m1, m2])
end

const Triple_Bond = begin
    m1 = Makie._mantle(Point3f(-4/5, 0.0, 1.0), Point3f(-4/5, 0.0, 1.0), _Triple_Radius, _Triple_Radius, _Bond_Quality)
    m2 = Makie._mantle(Point3f( 4/5, 0.0, 1.0), Point3f( 4/5, 0.0, 1.0), _Triple_Radius, _Triple_Radius, _Bond_Quality)
    m3 = Makie._mantle(Point3f( 0.0, 0.0, 0.0), Point3f( 0.0, 0.0, 1.0), _Triple_Radius, _Triple_Radius, _Bond_Quality)
    merge([m1, m2, m3])
end

const Aromatic_Bond = begin
    # 点線？
    m1 = Makie._mantle(Point3f(zeros(3)), Point3f((0,0,1)), _Bond_Radius, _Bond_Radius, _Bond_Quality)
    m2 = Makie._mantle(Point3f(zeros(3)), Point3f((0,0,1)), _Bond_Radius, _Bond_Radius, _Bond_Quality)
    merge([m1, m2])
end
