function testADDTable(a, b)
    c = {}
    c[1] = a[1] + b[1]
    c[2] = a[2] + b[2]
    c[3] = a[3] + b[3]
    c[4] = a[4] + b[4]
    return c
end

function testMAPTable(a, froimmin, frommax, tomin, tomax)
    c = {}
    c[1] = ((a[1] - froimmin) * ((tomax - tomin) / (frommax - froimmin))) + tomin
    c[2] = ((a[2] - froimmin) * ((tomax - tomin) / (frommax - froimmin))) + tomin
    c[3] = ((a[3] - froimmin) * ((tomax - tomin) / (frommax - froimmin))) + tomin
    c[4] = ((a[4] - froimmin) * ((tomax - tomin) / (frommax - froimmin))) + tomin
    return c
end


first_time = os.clock()

for i = 1, 10 do

    local a = { 1, 2, 3, 4 }
    local b = { -20.0, 100.0, 80.0, 200.0 }
    local c = addVec4(a, b)
    local c = mapVec4(212.0, 32.0, 100.0, 0.0, c) -- farenheit to celsius
    --print(c[1].." "..c[2].." "..c[3].." "..c[4])
end

last_time = os.clock()

print(last_time - first_time)

first_time = os.clock()

for i = 1, 10 do

    local a = { 1, 2, 3, 4 }
    local b = { -20.0, 100.0, 80.0, 200.0 }
    local c = testADDTable(a, b)
    local c = testMAPTable(c, 212.0, 32.0, 100.0, 0.0) -- farenheit to celsius
    --print(c[1].." "..c[2].." "..c[3].." "..c[4])
end

last_time = os.clock()

print(last_time - first_time)

print("enf of lua")