
"""
API for fetching consultation records and phenotype list
"""
def try_fix_json(raw_str):
    cleaned = re.split(r'</html>|<html', raw_str, flags=re.IGNORECASE)[0].strip()

    def fix_missing_commas(s):
        # 在結尾物件或陣列後若接字串（通常是 key）則補逗號
        s = re.sub(r'([}\]])\s*"(?=\w+"\s*:)', r'\1, "', s)
        # 一般物件/陣列之間直接連在一起也補逗號
        s = re.sub(r'([}\]])\s*(?=[{\[])', r'\1, ', s)
        return s

    def remove_trailing_commas(s):
        # 移除最後逗號前接閉合括號的情況
        s = re.sub(r',\s*([}\]])', r'\1', s)
        return s

    def balance_brackets(s):
        stack = []
        result = []
        quotes_open = False
        escape = False

        for i, char in enumerate(s):
            if char == '\\' and not escape:
                escape = True
                result.append(char)
                continue

            if char == '"' and not escape:
                quotes_open = not quotes_open

            if not quotes_open:
                if char in '{[':
                    stack.append(char)
                elif char in '}]':
                    if not stack:
                        continue
                    last = stack.pop()
                    if (last == '{' and char != '}') or (last == '[' and char != ']'):
                        stack.append(last)
                        continue

            escape = False
            result.append(char)

        while stack:
            open_br = stack.pop()
            result.append('}' if open_br == '{' else ']')

        return ''.join(result)

    # apply steps
    cleaned = fix_missing_commas(cleaned)
    cleaned = remove_trailing_commas(cleaned)
    cleaned = balance_brackets(cleaned)
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError as e:
        return f"Fail to fix the json：{e}"
    
def get_phenotype_list_from_api(data:dict):
    '''
    This API is designed for requesting the gene panel for a given medical record ID.
    Due to the incorrect design of this API, the results are returned in response.text instead of response.json. 
    Additionally, when multiple IDs are requested, the results are not separated into distinct dictionaries, 
    causing the final JSON to contain only the last sample.
    It is recommend to request this API by one ID.
    The following format is the expecting input 
    data={"ChartNo":"19946981"}
    '''

    url = "http://hisweb.hosp.ncku/hisservice/opd/nckuhisweb/aspx/DelegateExamServiceGate.aspx/GetPhenotypeList"
    
    payload = {
        "JasonInputValue": json.dumps(data)
    }
    response = requests.post(url, data=payload)
    text = response.text

    #aa=clean_broken_json(text)
    #return(aa)
    #return(text)
    
    json_part = text.replace('\r','').split("\n\n")[0]
    json_part = json_part.replace('\n','\\n')
    
    try:
        return(json.loads(json_part,strict=False))
    except json.JSONDecodeError as e:
        return(try_fix_json(json_part))
    
def get_consultation_record_from_api(data):
    '''
    This API is designed for requesting the consultation records for a given medical record ID.
    Multiple IDs are available.
    here is an expecting input format:
    data={"chartNo": "20303333",
          "tcode": "EMR-3-GC-002"}
    '''
    url = "https://apigw-i.apim.hosp.ncku.edu.tw/rd/prod-i/easyform/getdata"
    try:
        # Convert Python dictionary to JSON
        json_data = json.dumps(data)
        
        # Set appropriate headers
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json",
            "X-IBM-Client-Id": "9c03b0c83c562ffa22d1b4ff0e54d41d",
        }
        
        # Make the POST request
        response = requests.post(url, 
                               headers=headers,
                               data=json_data)
        
        # Check if the request was successful
        response.raise_for_status()
        
        # Return the JSON response
        return response.json()
    except requests.exceptions.HTTPError as errh:
        print(f"HTTP Error: {errh}")
    except requests.exceptions.ConnectionError as errc:
        print(f"Error connecting: {errc}")
    except requests.exceptions.Timeout as errt:
        print(f"Timeout Error: {errt}")
    except requests.exceptions.RequestException as err:
        print(f"Something went wrong: {err}")
    except json.JSONDecodeError:
        print("Failed to parse JSON response")
        return None

def fetch_medical_record(request,medical_record_number):
    '''
    data={"chartNo": medical_record_number,
        "tcode": "EMR-3-GC-002"}
    '''
    print(f'Fetch medical record for {medical_record_number}')
    genderOptions={'男':'M','女':'F'}
    if request.method == 'GET':
        consultation_record=get_consultation_record_from_api({"chartNo": medical_record_number,"tcode": "EMR-3-GC-002"})
        phenotype_list=get_phenotype_list_from_api({"ChartNo": medical_record_number})

        '''
        [{}] would be returned if nothing is fetched from get_phenotype_list_from_api 
        [] would be returned if nothing is fetched from get_consultation_record_from_api 
        '''
        if len(phenotype_list[0])==0 and len(consultation_record)==0:
            print(f'No record fetched for {medical_record_number}')
            result={}
        else:
            result={
                "gender": genderOptions[consultation_record[0].get('gender').replace(' ','')],
                "date_of_birth": consultation_record[0].get('date_of_birth'),
                "phenotype_list": [{
                    'date':phenotype_list[0].get('date'),
                    'phenotypes':phenotype_list[0].get('phenotypes')
                }] if len(phenotype_list[0])>0 else [],
                "consultation_record": [{
                    'consult':consultation_record[0].get('consult')
                }] if len(consultation_record)>0 else []
            }
        #print(result) ## for debugging
    else:
        result={}
    
    return JsonResponse(result)