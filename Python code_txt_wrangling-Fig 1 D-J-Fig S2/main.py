
import openpyxl


def write_data_to_xlsx(data):
    wb = openpyxl.Workbook()
    ws = wb.active
    for row in data:
        ws.append(row)
    wb.save("example.xlsx")
    print("success")


if __name__ == '__main__':
    data =[]
    with open("./New Text Document.txt", mode="r", encoding="utf-8") as file:
        for line in file:
            row = line.strip().split()
            if not row:
                continue
            if "val" in line:
                data.append([line.strip()])
            elif "Columns" in line:
                continue
            else:
                data.append(row)

    write_data_to_xlsx(data)
