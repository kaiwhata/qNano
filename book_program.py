from tkinter import *

class book:
    def __init__(self):
        title = ""




class GUI:
    def __init__(self):
        #def setupGUI():
        #global recordlist, name_field, author_field, ready_to_write, error_label

        self.canvas = Canvas(window, height = 0, width=400, bg='blue')
        self.canvas.create_rectangle(10, 20, 30, 10, fill='white')
        self.canvas.pack()


        name_label = Label(window, text='Enter name:')
        name_label.pack()
        self.name_field = Entry(window)
        self.name_field.pack()

        author_label = Label(window, text='Enter author:')
        author_label.pack()
        self.author_field = Entry(window)
        self.author_field.pack()

        button_label = Label(window, text='Press to validate:')
        button = Button(window, text='Submit', command=self.doSubmit)
        button_label1 = Label(window, text='Convert Record to csv')
        button1 = Button(window, text='write to csv', command=self.writetocsv)
        button_label.pack()
        button.pack()
        button_label1.pack()
        button1.pack()


        error_label = Label(window, text='')
        error_label.pack()

        #slider = Scale(window, from_=1, to=100, orient=HORIZONTAL, command=doSlide)
        #slider.pack()

        self.recordlist = [['bing', 'bill'], ['bam', 'grick'], ['boom', 'toffee']]

        self.ready_to_write = False
        window.mainloop()


    def doSubmit(self):
        #canvas.create_text(200, 200, text=name_field.get())
        #canvas.update()
        print(self.name_field.get())
        recordlist.append([self.name_field.get(),self.author_field.get() ])
        #author_field.delete(ALL) - command to clear field

        #data validation goes in here using ready_to_write variable

    def writetocsv(self):
        import csv
        if self.ready_to_write:
            ofile = open('database.txt', 'w') #open with write privelages
            writer = csv.writer(ofile, delimiter=',')
            for record in self.recordlist:
                print(record)
                writer.writerow(record)
            ofile.close()
        else:
            txt_var = StringVar()
            label_object = Label(window, textvariable=txt_var)
            label_object.pack()
            txt_var.set('Error, you need to Validate your data')
            #error_label.set('you havent validated your data') #Elf to check .set() command for labels




#def doSlide():
#    return

#setupGUI()
window = Tk()
window.title("Data Entry Screen")
GUI()
