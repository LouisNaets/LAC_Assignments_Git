

# Calculate and print partial safety factors for each load channel
load_channels = [
    {"name": "TbFA", "extreme": 215879.984375, "characteristic": 291437.97890625003, "design": 364297.47363281256},
    {"name": "TbSS", "extreme": 69760.18684895833, "characteristic": 94176.25224609375, "design": 117720.31530761719},
    {"name": "YbTilt", "extreme": 32585.416666666668, "characteristic": 43990.31250000001, "design": 54987.89062500001},
    {"name": "YbRoll", "extreme": 14167.382486979166, "characteristic": 19125.966357421876, "design": 23907.457946777344},
    {"name": "ShftTrs", "extreme": -12068.249674479166, "characteristic": 16292.137060546875, "design": 20365.171325683594},
    {"name": "OoPBRM", "extreme": -42829.711588541664, "characteristic": 57820.11064453125, "design": 72275.13830566406},
    {"name": "IPBRM", "extreme": 24160.418619791668, "characteristic": 32616.565136718753, "design": 40770.70642089844},
]

for channel in load_channels:
    SF1 = channel["design"] / channel["characteristic"]
    SF2 = channel["characteristic"] / channel["extreme"]
    print(f"Load channel: {channel['name']} [kNm]")
    print(f"Extreme value: {channel['extreme']}")
    print(f"Characteristic value: {channel['characteristic']}")
    print(f"Design value: {channel['design']}")
    print(f"Partial safety factor 1: {SF1}\n")
    print(f"Partial safety factor 2: {SF2}\n")
